#ifndef spiraltraj_gadgetron_H
#define spiraltraj_gadgetron_H

#include "vdspiral.h"

#include <ismrmrd/xml.h>
#include <gadgetron/log.h>
#include <gadgetron/Gadget.h>
#include <gadgetron/mri_core_utility.h>
#include "armadillo"
#include <gadgetron/mri_core_girf_correction.h>
#include <fstream>


namespace Gadgetron {
namespace Spiral {

    class spiraltraj_gadgetron {
    public:
        spiraltraj_gadgetron() = default;
        spiraltraj_gadgetron(const ISMRMRD::IsmrmrdHeader &h);

       std::pair<hoNDArray<float>, hoNDArray<float>>
        calculate_trajectories_and_weight(long ADCsamples);
        bool setDelayParams(const double& GDelay,const double& ADCShift);
        float getKmax() const;

    private:
        Core::optional<hoNDArray<std::complex<float>>> girf_kernel;
        float girf_sampling_time_us;
        long Tsamp_ns_;
        long Nints_;
        double gmax_;
        double smax_;
        double krmax_;
        double fov_;
        float TE_;
        //piv edit
        long Spiral_type;
        double Resolution_mm;
        double ADCShift_us;
        double GradDelay_us;
        double SpiralOS;
        vdspiral m_SpiralTraj;

        hoNDArray<floatd2> correct_gradients(const hoNDArray<floatd2> &gradients, float grad_samp_us,
                                             float girf_samp_us, const float *read_dir, const float *phase_dir,
                                             const float *slice_dir);


    };
}
}


#endif