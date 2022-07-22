#include "RemoveRefOSGadget.h"
#include "hoNDFFT.h"
#include "ismrmrd/xml.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron{

    RemoveRefOSGadget::RemoveRefOSGadget()
    {
    }

    RemoveRefOSGadget::~RemoveRefOSGadget()
    {
    }

    int RemoveRefOSGadget
        ::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
    {


            for(std::vector<IsmrmrdReconBit>::iterator it = m1->getObjectPtr()->rbit_.begin();
        it != m1->getObjectPtr()->rbit_.end(); ++it)
    {
        //Grab a reference to the buffer containing the imaging data
        //We are ignoring the reference data
        //IsmrmrdDataBuffered & dbuff2 = it->data_;
           
        IsmrmrdReconBit & rbit = *it;
        if (!rbit.ref_) continue;
    
        IsmrmrdDataBuffered & dbuff = (*rbit.ref_);
         
        std::vector<size_t> data_out_dims =dbuff.data_.dimensions();
        
       
    
    

        
        if ( !ifft_buf_.dimensions_equal(&data_out_dims) )
        {
            ifft_buf_.create(data_out_dims);
            ifft_res_.create(data_out_dims);
        }

        float ratioFOV = 2;

        data_out_dims[0] = (size_t)(data_out_dims[0]/ratioFOV);
        if ( !fft_buf_.dimensions_equal(&data_out_dims) )
        {
            fft_buf_.create(data_out_dims);
            fft_res_.create(data_out_dims);
        }

        size_t sRO = dbuff.data_.get_size(0); 

        hoNDFFT<float>::instance()->ifft1c(dbuff.data_, ifft_res_, ifft_buf_);

        vector_td<size_t, 4> crop_size{data_out_dims[0],data_out_dims[1],data_out_dims[2],data_out_dims[3]};
        vector_td<size_t, 4> crop_offset(0);
        crop_offset[0] = (size_t)( (sRO-data_out_dims[0]) / 2 );


        Gadgetron::crop(crop_offset,crop_size,ifft_res_,fft_buf_);
        hoNDFFT<float>::instance()->fft1c(fft_buf_,fft_res_);

      
        dbuff.data_=fft_res_;
        
        for (auto hdr: dbuff.headers_)
        {
            hdr.number_of_samples=hdr.number_of_samples/2;
        }

        
    }

      if (this->next()->putq(m1) == -1)
      {
        GERROR("RemoveRefOSGadget::process, passing data on to next gadget");
        return GADGET_FAIL;
      }

      return GADGET_OK;
    }


    GADGET_FACTORY_DECLARE(RemoveRefOSGadget)
}
