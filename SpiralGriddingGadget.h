#ifndef SpiralGriddingGadget_H
#define SpiralGriddingGadget_H

#include "Gadget.h"
#include "gadgetron_mricore_export.h"
#include "mri_core_data.h"


#include <gadgetron/hoNDFFT.h>
#include <gadgetron/hoNFFT.h>
#include <gadgetron/mri_core_utility.h>
#include <gadgetron/hoNDArray_utils.h>
#include <gadgetron/hoNDArray_elemwise.h>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE SpiralGriddingGadget : 
  public Gadget1<IsmrmrdReconData>
    {
    public:
      GADGET_DECLARE(SpiralGriddingGadget)
      SpiralGriddingGadget();


      GADGET_PROPERTY(use_calculated_csm, bool, "flag to use coil sens calculation in reconbit from bart", false);

      GADGET_PROPERTY(kernel_width,float,"Kernel width for NFFT", 5.5);
		  GADGET_PROPERTY(gridding_os,float,"Oversampling used in NFFT", 1.5);


    GADGET_PROPERTY(output_folder, std::string, "path of output folder for traj and pics", "/tmp/gadgetron/");
    GADGET_PROPERTY(perform_pics, bool,"Perform bart pics recon", false);
    GADGET_PROPERTY(pics_settings, std::string, "pics call", "pics -i 10 -l2 -r 1e-2");
    GADGET_PROPERTY(pics_scale_factor,float,"additional scaling of output image" , 1e3);

		
	//	GADGET_PROPERTY(image_series,int,"Image Series",1);




    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
      virtual int process_config(ACE_Message_Block* mb);
      long long image_counter_;
      ISMRMRD::MatrixSize matrix_size;
      ISMRMRD::IsmrmrdHeader header;


    };
}
#endif //SpiralGriddingGadget_H
