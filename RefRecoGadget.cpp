#include "RefRecoGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"
#include "CFL_IO.h"
namespace Gadgetron
{

    RefRecoGadget::RefRecoGadget()
        : image_counter_(0)
    {
    }

    int RefRecoGadget::process(GadgetContainerMessage<IsmrmrdReconData> *m1)
    {

        // Iterate over all the recon bits
        for (std::vector<IsmrmrdReconBit>::iterator it = m1->getObjectPtr()->rbit_.begin();
             it != m1->getObjectPtr()->rbit_.end(); ++it)
        {
            // Grab a reference to the buffer containing the imaging data
            // We are ignoring the reference data
            // IsmrmrdDataBuffered & dbuff2 = it->data_;

            IsmrmrdReconBit &rbit = *it;
            IsmrmrdDataBuffered &dbuff = (*rbit.ref_);
            if (!rbit.ref_)
                continue;

            std::array<uint16_t, 3> encode_matSz = rbit.data_.sampling_.encoded_matrix_;

            Gadgetron::uint64d4 ref_dim{encode_matSz[0], encode_matSz[1], encode_matSz[2], dbuff.data_.get_size(3)};
            auto buff_pad = hoNDArray<complex_float_t>(ref_dim[0], ref_dim[1], ref_dim[2], ref_dim[3]);
            auto buff_pad2 = hoNDArray<complex_float_t>(ref_dim[0], ref_dim[1], ref_dim[2], ref_dim[3]);

            Gadgetron::pad<complex_float_t, 4>(ref_dim, dbuff.data_, buff_pad2);

            Gadgetron::permute(buff_pad2,buff_pad,{1,0,2,3});

            // Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
            uint16_t E0 = buff_pad.get_size(0);
            uint16_t E1 = buff_pad.get_size(1);
            uint16_t E2 = buff_pad.get_size(2);
            uint16_t CHA = buff_pad.get_size(3);
            uint16_t N = buff_pad.get_size(4);
            uint16_t S = buff_pad.get_size(5);
            uint16_t LOC = buff_pad.get_size(6);

            Gadgetron::hoNDArray<complex_float_t> csm1;
            if (calc_csm.value())
            {
                
                // write calib data and do csm
                std::string calib_filename = output_folder.value() + "calibCFL";
                std::string sens_filename = output_folder.value() + "sensCFL";
                CFL_IO::hoNDArray2CFL<complex_float_t>(calib_filename, buff_pad);
                std::stringstream ss;
                ss.clear();
                ss << "bart ecalib " << ecalib_settings.value()<< " ";
                ss << calib_filename << " " << sens_filename << " ";
                std::string bart_cmd = ss.str();
                GDEBUG_STREAM("Executing : " << bart_cmd.c_str());
                std::system(bart_cmd.c_str());
                auto csm = CFL_IO::CFL2hoNDARRAY<complex_float_t>(sens_filename);

                csm1 = Gadgetron::hoNDArray<complex_float_t>(csm.get_size(0), csm.get_size(1), csm.get_size(2), csm.get_size(3));
                // overkill for option to flip csm. we can do it much nicer. But later
                uint32_t nVoxels = E0 * E1 * E2;
                for (uint16_t cha = 0; cha < CHA; cha++)
                {
                    Gadgetron::hoNDArray<complex_float_t>::iterator from_source = csm.begin() + (nVoxels * cha);
                    Gadgetron::hoNDArray<complex_float_t>::iterator from_target = csm1.begin() + (nVoxels * cha);
                    std::copy(from_source, from_source + nVoxels, from_target);
                }
                
                //replace reference data with coil sensitivities
                dbuff.data_=csm1;
            }
            
            if (export_reference_image.value())
            {
                // Create an image array message
                GadgetContainerMessage<IsmrmrdImageArray> *cm1 =
                    new GadgetContainerMessage<IsmrmrdImageArray>();

                // Grab references to the image array data and headers
                IsmrmrdImageArray &imarray = *cm1->getObjectPtr();

                // The image array data will be [E0,E1,E2,1,N,S,LOC] big
                // Will collapse across coils at the end
                std::vector<size_t> data_dims(7);
                data_dims[0] = E0;
                data_dims[1] = E1;
                data_dims[2] = E2;
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

                // We will not add any meta data
                // so skip the meta_ part

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
                            imarray.headers_(n, s, loc).matrix_size[0] = E0;
                            imarray.headers_(n, s, loc).matrix_size[1] = E1;
                            imarray.headers_(n, s, loc).matrix_size[2] = E2;
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

                            // Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
                            // Each chunk will be [E0,E1,E2,CHA] big
                            std::vector<size_t> chunk_dims(4);
                            chunk_dims[0] = E0;
                            chunk_dims[1] = E1;
                            chunk_dims[2] = E2;
                            chunk_dims[3] = CHA;
                            hoNDArray<std::complex<float>> chunk = hoNDArray<std::complex<float>>(chunk_dims, &buff_pad(0, 0, 0, 0, n, s, loc));

                            // Do the FFTs in place
                            hoNDFFT<float>::instance()->ifft3c(chunk);

                            // Square root of the sum of squares
                            // Each image will be [E0,E1,E2,1] big
                            std::vector<size_t> img_dims(3);
                            img_dims[0] = E0;
                            img_dims[1] = E1;
                            img_dims[2] = E2;
                            hoNDArray<std::complex<float>> output = hoNDArray<std::complex<float>>(img_dims, &imarray.data_(0, 0, 0, 0, n, s, loc));
                            // Zero out the output
                            clear(output);

                            if (calc_csm.value())
                            { // adpative coil combine
                                multiplyConj(chunk, csm1, chunk);
                                // Add up
                                for (size_t c = 0; c < CHA; c++)
                                {
                                    output += hoNDArray<std::complex<float>>(img_dims, &chunk(0, 0, 0, c));
                                }
                            }
                            else // sos
                            {
                                // Compute d* d in place
                                multiplyConj(chunk, chunk, chunk);
                                // Add up
                                for (size_t c = 0; c < CHA; c++)
                                {
                                    output += hoNDArray<std::complex<float>>(img_dims, &chunk(0, 0, 0, c));
                                }
                                // Take the square root in place
                                sqrt_inplace(&output);
                            }
                        }
                    }
                }

                // Pass the image array down the chain
                if (this->next()->putq(cm1) < 0)
                {
                    m1->release();
                    return GADGET_FAIL;
                }
            } // if (export_reference_image)
        }     // reconbit iteration

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("Put IsmrmrdReconData to Q failed ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    } // process()

    GADGET_FACTORY_DECLARE(RefRecoGadget)
}
