#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include "mri_core_data.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

    class EXPORTGADGETSMRICORE RemoveRefOSGadget :
        public Gadget1<IsmrmrdReconData>
    {
    public:
        GADGET_DECLARE(RemoveRefOSGadget);

        RemoveRefOSGadget();
        virtual ~RemoveRefOSGadget();

    protected:

        virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);

        hoNDArray< std::complex<float> > fft_res_;
        hoNDArray< std::complex<float> > ifft_res_;

        hoNDArray< std::complex<float> > fft_buf_;
        hoNDArray< std::complex<float> > ifft_buf_;

        int   encodeNx_;
        float encodeFOV_;
        int   reconNx_;
        float reconFOV_;

    };
}
