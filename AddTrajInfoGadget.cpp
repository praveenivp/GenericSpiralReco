#include "AddTrajInfoGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"

namespace Gadgetron
{

    AddTrajInfoGadget::AddTrajInfoGadget()
        : image_counter_(0)
    {
    }

    int AddTrajInfoGadget::process_config(ACE_Message_Block *mb)
    {

        try
        {
            deserialize(mb->rd_ptr(), MeasHeader);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!MeasHeader.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        SpiralTraj = Spiral::spiraltraj_gadgetron(MeasHeader);
        SpiralTraj.setDelayParams(grad_delay.value(),ADC_shift.value());

        return GADGET_OK;
    }

    int AddTrajInfoGadget::process(GadgetContainerMessage<IsmrmrdReconData> *m1)
    {

        // Iterate over all the recon bits
        for (std::vector<IsmrmrdReconBit>::iterator it = m1->getObjectPtr()->rbit_.begin();
             it != m1->getObjectPtr()->rbit_.end(); ++it)
        {
            // Grab a reference to the buffer containing the imaging data
            // We are ignoring the reference data
            // IsmrmrdDataBuffered & dbuff2 = it->data_;

            IsmrmrdReconBit &rbit = *it;
            IsmrmrdDataBuffered &dbuff = it->data_;

            
            // trajectories are -0.5 to 0.5 scaled  with dim (Grad_Axis=2,ADC_points, Interleaves/E1)
            auto traj = SpiralTraj.calculate_trajectories_and_weight(dbuff.headers_[0].number_of_samples);
            GDEBUG("done calculating trajectories\n");

           
            dbuff.trajectory_ = traj.first;
            dbuff.density_ = traj.second;

            if (calc_3D_traj.value())
            {
                // trajectories are -0.5 to 0.5 scaled  with dim (Grad_Axis=3,ADC_points, Interleaves/E1,Partitions/E2)
                CalcTrajectory3D(dbuff);
            }
            
         RemoveUnacquiredData(dbuff);

            PerformFOVShift(dbuff.data_, dbuff.headers_[0], dbuff.trajectory_.value());



        } // reconbit iteration

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("Put IsmrmrdReconData to Q failed ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    } // process()

    bool AddTrajInfoGadget::PerformFOVShift(Gadgetron::hoNDArray<complex_float_t> &ksp_data, ISMRMRD::AcquisitionHeader &acq_hdr, Gadgetron::hoNDArray<float> &kTraj)
    {
          using namespace Gadgetron::Indexing;
        // ksp_data 2 dimensional COLx CHA 
        // kTRaj 2 dimensional COLxLIN
        const float pi = std::acos(-1);
        float kmax = (2.f * pi) * SpiralTraj.getKmax(); // 1/m
        auto sz = ksp_data.dimensions();
        auto sz1 = kTraj.dimensions();

        auto slcPos = acq_hdr.position; // mm
                                        // if (acq_hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_REPETITION))
        GDEBUG_STREAM("RPOS " << slcPos[0] << " PhPos " << slcPos[1] << " slcPos " << slcPos[2]);

      


        size_t E0 = ksp_data.get_size(0);
        size_t E1 = ksp_data.get_size(1);
        size_t E2 = ksp_data.get_size(2);
        size_t CHA = ksp_data.get_size(3);
        size_t N = ksp_data.get_size(4);
        size_t S = ksp_data.get_size(5);
        size_t LOC = ksp_data.get_size(6);
        


        //calculate the signal phase modulation because of off-center shift
        const std::complex<float> imagi(0, 1); // 0+1i
        Gadgetron::hoNDArray<complex_float_t> B0mod(kTraj.get_size(1), kTraj.get_size(2));
        size_t traj_dim=kTraj.get_size(0);
        size_t ELEM = E0 * E1 *traj_dim;
        for (uint32_t i = 0; i < ELEM; i += traj_dim)
        {
            B0mod[i / traj_dim] = (imagi * ((kTraj[i + 1] * slcPos[0] * 1e-3f + kTraj[i] * slcPos[1] * -1e-3f) * (2.f * kmax)));
        }
        // std::transform(test.begin(),test.end(),B0mod.begin(),[kmax,slcPos,pi,imagi](auto x){return (imagi*((2.f*pi)*(x[0]*slcPos[0]*1e-3f+x[1]*slcPos[1]*-1e-3f)*(2.f*kmax)));});

 
        // apply that off-center shift to the data;
        for (uint16_t loc = 0; loc < LOC; loc++)
        {
            for (uint16_t s = 0; s < S; s++)
            {
                for (uint16_t n = 0; n < N; n++)
                {
                    for (uint16_t cha = 0; cha < CHA; cha++)
                    {
                        for (uint16_t e2 = 0; e2 < E2; e2++)
                        {
                            complex_float_t *startPtr = &ksp_data(0, 0, e2, cha, n, s, loc);
                            for (auto elem = 0; elem < E0 * E1; ++elem)
                            {
                                startPtr[elem] = startPtr[elem] * std::exp(B0mod[elem]);
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    bool AddTrajInfoGadget::RemoveUnacquiredData(IsmrmrdDataBuffered &dbuff)
    {

        if (!MeasHeader.encoding[0].parallelImaging.has_value())
            return false;

        auto R_E1 = MeasHeader.encoding[0].parallelImaging.value().accelerationFactor.kspace_encoding_step_1;
        auto R_E2 = MeasHeader.encoding[0].parallelImaging.value().accelerationFactor.kspace_encoding_step_2;

        if(R_E1*R_E2 ==1)
        return false;
        
        GDEBUG_STREAM("Trying to remove unacquired datapoints: R_E1= "<<R_E1 <<" ,R_E2= "<<R_E2<<std::endl);

        auto &h = dbuff.headers_;
        auto &traj = dbuff.trajectory_.value();
        auto &dcw = dbuff.density_.value();
        auto &data = dbuff.data_;
        

        size_t E0 = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);
        size_t N = data.get_size(4);
        size_t S = data.get_size(5);
        size_t LOC = data.get_size(6);

        auto new_headers=hoNDArray<ISMRMRD::AcquisitionHeader>(E1/R_E1,E2/R_E2,N,S,LOC);
        auto new_data=hoNDArray<complex_float_t>(E0,E1/R_E1,E2/R_E2,CHA,N,S,LOC);
        auto new_traj=hoNDArray<float>(3,E0,E1/R_E1,E2/R_E2);
        auto new_dcw=hoNDArray<float>(E0,E1/R_E1,E2/R_E2);

        std::vector<uint16_t> E1_to_copy;
        std::vector<uint16_t> E2_to_copy;
        for (uint16_t loc = 0; loc < LOC; loc++)
        {
            for (uint16_t s = 0; s < S; s++)
            {
                for (uint16_t n = 0; n < N; n++)
                {
                    for (uint16_t e2 = 0; e2 < E2; e2++)
                    {
                        for (uint16_t e1 = 0; e1 < E1; e1++)
                        {

                            auto & nsamp=h(e1,e2,n,s,loc).number_of_samples;
                            if(nsamp>0)
                            {
                                E1_to_copy.push_back(h(e1,e2,n,s,loc).idx.kspace_encode_step_1);
                                E2_to_copy.push_back(h(e1,e2,n,s,loc).idx.kspace_encode_step_2);
                                //CAIPI _SHIFT will fail
                            new_headers(e1/R_E1,e2/R_E2,n,s,loc)=h(e1,e2,n,s,loc);

                            for (uint16_t cha=0; cha<CHA;cha++)
                            {
                                complex_float_t* FromPtr= &data(0,e1,e2,cha,n,s,loc);
                                complex_float_t* ToPtr= &new_data(0,e1/R_E1,e2/R_E2,cha,n,s,loc);
                                std::copy(FromPtr,FromPtr+E0,ToPtr);
                            }
                            }

                        }
                    }
                        
                }
            }
        }

                    auto traj_dim=traj.get_size(0);
                    float* ToPtr_traj= &new_traj[0];
                    float* ToPtr_dcw= &new_dcw[0];
                    for (uint16_t idx = 0; idx <E1_to_copy.size() ; idx++)
                    {
                      
                                //traj
                                float* FromPtr_traj= &traj(0,0,E1_to_copy[idx],E2_to_copy[idx]);
                                
                                std::copy(FromPtr_traj,FromPtr_traj+(E0*traj_dim),ToPtr_traj);

                                   //dcw
                                float* FromPtr= &dcw(0,E1_to_copy[idx],E2_to_copy[idx]);
                                std::copy(FromPtr,FromPtr+E0,ToPtr_dcw);
                                ToPtr_traj+=(E0*traj_dim);
                                ToPtr_dcw+=(E0);

                    }
                            

        dbuff.data_=new_data;
        dbuff.density_=new_dcw;
        dbuff.headers_=new_headers;
        dbuff.trajectory_=new_traj;

        GDEBUG("Sucssfully modified recon buff\n");
        return true;
    }

    bool AddTrajInfoGadget::CalcTrajectory3D(IsmrmrdDataBuffered &dbuff)
    {
        // make the first dim 2 to 3
        auto & kTraj2D=dbuff.trajectory_.value();
        auto & dcw2D=dbuff.density_.value();
        auto mat_sz = MeasHeader.encoding[0].encodedSpace.matrixSize;
        auto ktraj3D = Gadgetron::hoNDArray<float>(3, kTraj2D.get_size(1), kTraj2D.get_size(2), mat_sz.z);
        auto dcw3D = Gadgetron::hoNDArray<float>(dcw2D.get_size(0), dcw2D.get_size(1), mat_sz.z);

        auto nPts_2D = kTraj2D.get_size(1) * kTraj2D.get_size(2);
        float kz = 0.0f;
        for (size_t z = 0; z < mat_sz.z; z++)
        {
            auto startPtr = &ktraj3D(0, 0, 0, z);
            if (mat_sz.z > 1)
                kz = static_cast<float>(z) / static_cast<float>(mat_sz.z) - 0.5f;

            for (size_t elem = 0; elem < nPts_2D; elem++)
            {
                startPtr[elem * 3] = kTraj2D[elem * 2];
                startPtr[elem * 3 + 1] = kTraj2D[elem * 2 + 1];
                startPtr[elem * 3 + 2] = kz;
            }

            auto startPtr_dcw=&dcw3D(0,0,z);
            std::copy(&dcw2D[0],&dcw2D[0]+dcw2D.get_number_of_elements(),startPtr_dcw);

        }
        dbuff.trajectory_=ktraj3D;
        dbuff.density_=dcw3D;


        return true;
    }
    GADGET_FACTORY_DECLARE(AddTrajInfoGadget)
}
