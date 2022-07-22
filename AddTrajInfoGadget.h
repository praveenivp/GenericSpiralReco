#ifndef AddTrajInfoGadget_H
#define AddTrajInfoGadget_H

#include "Gadget.h"
#include "gadgetron_mricore_export.h"

#include "mri_core_data.h"

#include "spiraltraj/spiraltraj_gadgetron.h"

namespace Gadgetron{

  class EXPORTGADGETSMRICORE AddTrajInfoGadget : 
  public Gadget1<IsmrmrdReconData>
    {
    public:
      GADGET_DECLARE(AddTrajInfoGadget)
      AddTrajInfoGadget();


      GADGET_PROPERTY(calc_csm, bool, "flag to enable bart coil sens calculation", true);
      GADGET_PROPERTY(export_reference_image, bool, "export coil combined image from reference data", true);

      GADGET_PROPERTY(grad_delay, float, "additional gradient delay in us", 0);
      GADGET_PROPERTY(calc_3D_traj, bool, "flag to calculate 3D kspace location", true);



    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
      virtual int process_config(ACE_Message_Block* mb);
      bool PerformFOVShift(Gadgetron::hoNDArray<complex_float_t> &ksp_data, ISMRMRD::AcquisitionHeader &acq_hdr, Gadgetron::hoNDArray<float> &kTraj);
      bool RemoveUnacquiredData(IsmrmrdDataBuffered &dbuff);
      bool CalcTrajectory3D(IsmrmrdDataBuffered &dbuff);
      long long image_counter_;
      ISMRMRD::IsmrmrdHeader MeasHeader;
    
    private:
      Spiral::spiraltraj_gadgetron SpiralTraj;


    };
}
#endif //AddTrajInfoGadget_H
