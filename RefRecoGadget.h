#ifndef RefRecoGadget_H
#define RefRecoGadget_H

#include "Gadget.h"
#include "gadgetron_mricore_export.h"

#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETSMRICORE RefRecoGadget : 
  public Gadget1<IsmrmrdReconData>
    {
    public:
      GADGET_DECLARE(RefRecoGadget)
      RefRecoGadget();


      GADGET_PROPERTY(calc_csm, bool, "flag to enable bart coil sens calculation", true);
      GADGET_PROPERTY(export_reference_image, bool, "export coil combined image from reference data", true);

      GADGET_PROPERTY(ecalib_settings, std::string, "ecalib call", "ecalib -m2 -k4 -r24");
      GADGET_PROPERTY(output_folder, std::string, "path of output folder for calib data and coil sens cfl files", "/tmp/gadgetron/");


    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
      long long image_counter_;

    };
}
#endif //RefRecoGadget_H
