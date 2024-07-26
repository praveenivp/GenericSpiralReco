#include "SpiralGriddingGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "CFL_IO2.h"

#include <gadgetron/mri_core_coil_map_estimation.h>

namespace Gadgetron
{

    SpiralGriddingGadget::SpiralGriddingGadget()
        : image_counter_(0)
    {
    }

    int SpiralGriddingGadget::process_config(ACE_Message_Block *mb)
    {

        try
        {
            deserialize(mb->rd_ptr(), header);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!header.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int SpiralGriddingGadget::process(GadgetContainerMessage<IsmrmrdReconData> *m1)
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

            auto &data = dbuff.data_;
            size_t E0 = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t LOC = data.get_size(6);

            matrix_size = header.encoding.front().encodedSpace.matrixSize;
            hoNFFT_plan<float, 2> nufft = hoNFFT_plan(vector_td<size_t, 2>(matrix_size.x, matrix_size.y), gridding_os.value(), kernel_width.value());
            // NFFT_prep_mode mode=NFFT_prep_mode::NC2C;

            if (!perform_pics.value())
            {

                auto &traj = dbuff.trajectory_.value();
                Gadgetron::hoNDArray<floatd2> trajectory(traj.get_size(1) * traj.get_size(2));
                for (uint32_t i = 0; i < E0 * E1 * traj.get_size(0); i += traj.get_size(0))
                {
                    trajectory[i / traj.get_size(0)] = vector_td<float, 2>(traj[i], traj[i + 1]);
                }
                nufft.preprocess(trajectory, NFFT_prep_mode::NC2C);
            }

            auto &dcw = dbuff.density_.value();
            dcw.reshape(-1);

            // Create an image array message
            GadgetContainerMessage<IsmrmrdImageArray> *cm1 =
                new GadgetContainerMessage<IsmrmrdImageArray>();

            // Grab references to the image array data and headers
            IsmrmrdImageArray &imarray = *cm1->getObjectPtr();

            // The image array data will be [E0,E1,E2,1,N,S,LOC] big
            // Will collapse across coils at the end
            std::vector<size_t> data_dims(7);
            data_dims[0] = matrix_size.x;
            data_dims[1] = matrix_size.y;
            data_dims[2] = matrix_size.z;
            data_dims[3] = 1;
            data_dims[4] = N;
            data_dims[5] = S;
            data_dims[6] = LOC;
            imarray.data_.create(data_dims);

            // ImageHeaders will be [N, S, LOC]
            std::vector<size_t> header_dims(3);
            header_dims[0] = N;
            header_dims[1] = S;
            header_dims[2] = LOC;
            imarray.headers_.create(header_dims);

            std::vector<size_t> img_dims(3);
            img_dims[0] = matrix_size.x;
            img_dims[1] = matrix_size.y;
            img_dims[2] = matrix_size.z;

            // We will not add any meta data
            // so skip the meta_ part

            auto &headers = dbuff.headers_;
            // Loop over S and N and LOC
            for (uint16_t loc = 0; loc < LOC; loc++)
            {
                for (uint16_t s = 0; s < S; s++)
                {
                    for (uint16_t n = 0; n < N; n++)
                    {



                        // Set some information into the image header
                        // Use the middle acquisition header for some info
                        //[E1, E2, N, S, LOC]
                        ISMRMRD::AcquisitionHeader &acqhdr = dbuff.headers_(0,
                                                                            0,
                                                                            n, s, loc);
                        imarray.headers_(n, s, loc).measurement_uid = acqhdr.measurement_uid;
                        imarray.headers_(n, s, loc).matrix_size[0] = matrix_size.x;
                        imarray.headers_(n, s, loc).matrix_size[1] = matrix_size.y;
                        imarray.headers_(n, s, loc).matrix_size[2] = matrix_size.z;
                        imarray.headers_(n, s, loc).field_of_view[0] = dbuff.sampling_.encoded_FOV_[0];
                        imarray.headers_(n, s, loc).field_of_view[1] = dbuff.sampling_.encoded_FOV_[1];
                        imarray.headers_(n, s, loc).field_of_view[2] = dbuff.sampling_.encoded_FOV_[2];
                        imarray.headers_(n, s, loc).channels = 1;
                        imarray.headers_(n, s, loc).average = acqhdr.idx.average;
                        imarray.headers_(n, s, loc).slice = acqhdr.idx.slice;
                        imarray.headers_(n, s, loc).contrast = acqhdr.idx.contrast;
                        imarray.headers_(n, s, loc).phase = acqhdr.idx.phase;
                        imarray.headers_(n, s, loc).repetition = acqhdr.idx.repetition;
                        imarray.headers_(n, s, loc).set = acqhdr.idx.set;
                        imarray.headers_(n, s, loc).acquisition_time_stamp = acqhdr.acquisition_time_stamp;
                        imarray.headers_(n, s, loc).position[0] = acqhdr.position[0];
                        imarray.headers_(n, s, loc).position[1] = acqhdr.position[1];
                        imarray.headers_(n, s, loc).position[2] = acqhdr.position[2];
                        imarray.headers_(n, s, loc).read_dir[0] = acqhdr.read_dir[0];
                        imarray.headers_(n, s, loc).read_dir[1] = acqhdr.read_dir[1];
                        imarray.headers_(n, s, loc).read_dir[2] = acqhdr.read_dir[2];
                        imarray.headers_(n, s, loc).phase_dir[0] = acqhdr.phase_dir[0];
                        imarray.headers_(n, s, loc).phase_dir[1] = acqhdr.phase_dir[1];
                        imarray.headers_(n, s, loc).phase_dir[2] = acqhdr.phase_dir[2];
                        imarray.headers_(n, s, loc).slice_dir[0] = acqhdr.slice_dir[0];
                        imarray.headers_(n, s, loc).slice_dir[1] = acqhdr.slice_dir[1];
                        imarray.headers_(n, s, loc).slice_dir[2] = acqhdr.slice_dir[2];
                        imarray.headers_(n, s, loc).patient_table_position[0] = acqhdr.patient_table_position[0];
                        imarray.headers_(n, s, loc).patient_table_position[1] = acqhdr.patient_table_position[1];
                        imarray.headers_(n, s, loc).patient_table_position[2] = acqhdr.patient_table_position[2];
                        imarray.headers_(n, s, loc).data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                        imarray.headers_(n, s, loc).image_index = ++image_counter_;
                        imarray.headers_(n, s, loc).image_series_index = 1;

                        if (perform_pics.value()) // do pics
                        {
                            std::string traj_filename = output_folder.value() + "Traj";
                            std::string DCF_filename = output_folder.value() + "DCF";
                            std::string data_filename = output_folder.value() + "coil_data";

                            // Export Trajectory
                            auto &traj_unscaled = dbuff.trajectory_.value();
                            hoNDArray<complex_float_t> traj_scaled = hoNDArray<complex_float_t>(traj_unscaled.get_size(0), traj_unscaled.get_size(1), traj_unscaled.get_size(2), traj_unscaled.get_size(3));

                            auto mat_sz = header.encoding[0].reconSpace.matrixSize;
                            for (uint32_t i = 0; i < E0 * E1 * E2 * traj_scaled.get_size(0); i += 3)
                            {
                                traj_scaled[i].real(traj_unscaled[i] * static_cast<float>(mat_sz.x));
                                traj_scaled[i + 1].real(traj_unscaled[i + 1] * static_cast<float>(mat_sz.y));
                                traj_scaled[i + 2].real(traj_unscaled[i + 2] * static_cast<float>(mat_sz.z));
                            }

                            traj_scaled.reshape(3, E0, E1 * E2);
                            CFL_IO2::hoNDArray2CFL2<complex_float_t>(traj_filename, traj_scaled);

                            // export density compensation
                            auto &dcf = dbuff.density_.value();
                            hoNDArray<complex_float_t> h0ND_dcf = hoNDArray<complex_float_t>(1, E0, E1 * E2);
                            for (uint32_t i = 0; i < E0 * E1 * E2; i++)
                            {
                                h0ND_dcf[i].real(std::sqrt(dcf[i]));
                            }
                            CFL_IO2::hoNDArray2CFL2<complex_float_t>(DCF_filename, h0ND_dcf);

                            // export Coildata

                            std::vector<size_t> chunk_dims2(5);
                            chunk_dims2[0] = 1;
                            chunk_dims2[1] = dbuff.data_.get_size(0);
                            chunk_dims2[2] = dbuff.data_.get_size(1) * dbuff.data_.get_size(2);
                            chunk_dims2[3] = CHA;
                            chunk_dims2[4] = 1;

                            hoNDArray<std::complex<float>> chunk2 = hoNDArray<std::complex<float>>(chunk_dims2, &dbuff.data_(0, 0, 0, 0, n, s, loc));

                            CFL_IO2::hoNDArray2CFL2<complex_float_t>(data_filename, chunk2);

                            // do bart pics

                            std::string sens_filename = output_folder.value() + "sensCFL"; // from RefRecoGadget
                            std::string img_filename = output_folder.value() + "pics_reco";

                            // bart pics -i 10 -l2 -r 1e-2 -p DCF  -t Traj coil_data sensCFL out1
                            std::stringstream ss;
                            ss.clear();
                            ss << "bart " << pics_settings.value() << " ";
                            ss << "-p"
                               << " " << DCF_filename << " ";
                            ss << "-t"
                               << " " << traj_filename << " ";
                            ss << data_filename << " " << sens_filename << " " << img_filename;
                            std::string bart_cmd = ss.str();
                            GDEBUG_STREAM("Executing : " << bart_cmd.c_str());
                            std::system(bart_cmd.c_str());
                            auto pics_im = CFL_IO2::CFL2hoNDARRAY2<std::complex<float>>(img_filename);

                            hoNDArray<std::complex<float>> output = hoNDArray<std::complex<float>>(img_dims, &imarray.data_(0, 0, 0, 0, n, s, loc));
                            uint32_t nVoxels = (img_dims[2] > 0) ? img_dims[0] * img_dims[1] * img_dims[2] : img_dims[0] * img_dims[1];
                            std::copy(pics_im.begin(), pics_im.begin() + nVoxels, output.begin());

                            output*=pics_scale_factor.value(); //auto scaling gadget sucks!

                        } 
                        else // do gridding and coil combination
                        {

                            // Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
                            // Each chunk will be [E0,E1,E2,CHA] big
                            std::vector<size_t> chunk_dims(4);
                            chunk_dims[0] = dbuff.data_.get_size(0);
                            chunk_dims[1] = dbuff.data_.get_size(1);
                            chunk_dims[2] = dbuff.data_.get_size(2);
                            chunk_dims[3] = CHA;

                            hoNDArray<std::complex<float>> chunk = hoNDArray<std::complex<float>>(chunk_dims, &dbuff.data_(0, 0, 0, 0, n, s, loc));

                            chunk.reshape(E0 * E1, E2 * CHA);
                            hoNDArray<std::complex<float>> coil_images(matrix_size.x, matrix_size.y, E2 * CHA);
                            nufft.compute(chunk, coil_images, &dcw, NFFT_comp_mode::BACKWARDS_NC2C);

                            coil_images.reshape(matrix_size.x, matrix_size.y, E2, CHA);

                            if (E2 > 1) // do 1D-fft along E2 if 3D acquistion
                            {
                                coil_images = permute(coil_images, {2, 0, 1, 3});
                                hoNDFFT<float>::instance()->ifft1c(coil_images);
                                coil_images = permute(coil_images, {1, 2, 0, 3});
                            }

                            // Square root of the sum of squares
                            // Each image will be [E0,E1,E2,1] big
                            hoNDArray<std::complex<float>> output = hoNDArray<std::complex<float>>(img_dims, &imarray.data_(0, 0, 0, 0, n, s, loc));
                            // Zero out the output
                            clear(output);

                            if (use_calculated_csm.value() && it->ref_.has_value()) // use coil sens from RefRecoGadget()
                            {                                                       // adpative coil combine
                                hoNDArray<complex_float_t> &csm = it->ref_.value().data_;
                                multiplyConj(coil_images, csm, coil_images);
                                // Add up
                                for (size_t c = 0; c < CHA; c++)
                                {
                                    output += hoNDArray<std::complex<float>>(img_dims, &coil_images(0, 0, 0, c));
                                }
                            }
                            else // sum of squares
                            {
                                // // Compute d* d in place
                                multiplyConj(coil_images, coil_images, coil_images);
                                // Add up
                                for (size_t c = 0; c < CHA; c++)
                                {
                                    output += hoNDArray<std::complex<float>>(img_dims, &coil_images(0, 0, 0, c));
                                }
                                // // Take the square root in place
                                sqrt_inplace(&output);
                            } // sum of squares
                        }     // do pics else
                    }
                }
            }

            // Pass the image array down the chain
            if (this->next()->putq(cm1) < 0)
            {
                m1->release();
                return GADGET_FAIL;
            }

            // use and acquistion header and do undersampling in trajectory?

        } // reconbit iteration

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("Put IsmrmrdReconData to Q failed ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    } // process()

    GADGET_FACTORY_DECLARE(SpiralGriddingGadget)
}
