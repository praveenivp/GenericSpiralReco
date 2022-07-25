#include "spiraltraj_gadgetron.h"

#include <gadgetron/hoNDArray_utils.h>
#include <gadgetron/hoNDArray_fileio.h>
#include <gadgetron/vector_td_utilities.h>


namespace Gadgetron{
    namespace Spiral{
        spiraltraj_gadgetron::spiraltraj_gadgetron(const ISMRMRD::IsmrmrdHeader &h)
        {
                        ISMRMRD::TrajectoryDescription traj_desc;

            if (h.encoding[0].trajectoryDescription) {
                traj_desc = *h.encoding[0].trajectoryDescription;
            } else {
                throw std::runtime_error("Trajectory description missing");
            }

            if (traj_desc.identifier != "PESpiral") {
                throw std::runtime_error("Expected trajectory description identifier 'PESpiral', not found.");
            }


            try {
                auto userparam_long = to_map(traj_desc.userParameterLong);
                auto userparam_double = to_map(traj_desc.userParameterDouble);
                //piv edit
                Spiral_type=userparam_long.at("SpiralType");
                Resolution_mm=userparam_double.at("Resolution_mm");
                ADCShift_us=0;
                Tsamp_ns_ = userparam_long.at("DwellTime_ns");
                Nints_ = userparam_long.at("interleaves");
                GradDelay_us=3.85;
                SpiralOS=userparam_double.at("SpiralOS");

                gmax_ = userparam_double.at("MaxGradient_mT_per_m");
                smax_ = userparam_double.at("Slewmax_mT_m_ms");
                krmax_ =  0.5/(Resolution_mm*1e-3);//1/m
                fov_ = h.encoding[0].reconSpace.fieldOfView_mm.x; //userparam_double.at("FOVCoeff_1_cm");
            } catch (std::out_of_range exception) {
                std::string s = "Missing user parameters: " + std::string(exception.what());
                throw std::runtime_error(s);

            }
            TE_ = h.sequenceParameters->TE->at(0);

            GDEBUG("GradDelay:                    %f\n", GradDelay_us);
            GDEBUG("ADCshift:                    %f\n", ADCShift_us);
            GDEBUG("smax:                    %f\n", smax_);
            GDEBUG("gmax:                    %f\n", gmax_);
            GDEBUG("Tsamp_ns:                %d\n", Tsamp_ns_);
            GDEBUG("Nints:                   %d\n", Nints_);
            GDEBUG("fov:                     %f\n", fov_);
            GDEBUG("krmax (1/m):                   %f\n", krmax_);
            GDEBUG("Resolution_mm:                   %f\n",Resolution_mm);
            GDEBUG("SpiralType:                   %d\n", Spiral_type);
            GDEBUG("SpiralOS:                   %f\n", SpiralOS);
            GDEBUG("GIRF kernel:             %d\n", bool(this->girf_kernel));

        }

        bool spiraltraj_gadgetron::setDelayParams(const double& GDelay,const double& ADCShift)
        {
            this->GradDelay_us+=GDelay;
            this->ADCShift_us=ADCShift;
            GDEBUG("new GradDelay(us):                    %f\n", this->GradDelay_us);
            GDEBUG("new ADCshift(us):                    %f\n", this->ADCShift_us);
            return true;
        }

        float spiraltraj_gadgetron::getKmax() const
        {
            return krmax_; //1/m
        }

                std::pair<hoNDArray<float>, hoNDArray<float>>
        spiraltraj_gadgetron::calculate_trajectories_and_weight(long ADCsamples ) {            	


            int nfov = 4;
		    std::vector<double> v_fov({fov_,fov_,fov_,fov_});
            if(SpiralOS>1.0) 
            {
                v_fov[0]*=SpiralOS;
                v_fov[1]*=SpiralOS;
                GDEBUG("center FOV over sampled:                    %f\n", v_fov[1]);
            }

            std::vector<double> v_radius({0,0.18,0.25,1});
		
            double dGradMaxAmpl = gmax_;
            double gammabar = 42.5766; //kHz/mT
            vdspiral::eSpiralType typ = vdspiral::eSpiralType(Spiral_type);
            
           
            double dMinRiseTime=1000./smax_;
            m_SpiralTraj.prep(Nints_,Resolution_mm,v_fov,v_radius,dGradMaxAmpl,dMinRiseTime,typ,gammabar,GRAD_RASTER_TIME);
       
            //long ADCsamples=(long)acq_header.number_of_samples;
            std::vector<float> kx((Nints_*ADCsamples),0.);
            std::vector<float> ky((Nints_*ADCsamples),0.);;
            std::vector<float> dcf((Nints_*ADCsamples),0.);;
            bool suceeded=m_SpiralTraj.calcTrajectory(kx,ky,dcf,ADCsamples,static_cast<uint16_t>(fov_/Resolution_mm),ADCShift_us,GradDelay_us);



            GDEBUG("ADCsamples:             %d\n", ADCsamples);

            // //write trajectories for debugging
            // std::ofstream myfile2;
            // myfile2.open ("/tmp/gadgetron/Trajectories.txt",std::ofstream::trunc);
            // myfile2 <<"kx [1/m]          "<<" ky[1/m]         "<<"DCF        "<<'\n';
            // for(auto i=0;i<kx.size();++i){
            //     myfile2 <<kx[i] <<"  "<<ky[i]<<"  "<<dcf[i]<<'\n';
            //     if((i+1)%ADCsamples==0) myfile2 <<"\n\n";
            // }
            // myfile2<<std::flush;
            // myfile2.close();

            // scale and pack trajectory 
            float fKmax=krmax_; //1/m
            hoNDArray<float> trajectories2(ADCsamples,Nints_,2);
            size_t nsamples = Nints_*ADCsamples;


           //std::transform(kx.begin(),kx.end(),ky.begin(),trajectories.begin(),[fKmax](auto x, auto y){return floatd2(y,x)/(2*fKmax);});
           //inplace scaling
           std::transform(kx.begin(),kx.end(),trajectories2.begin(),[fKmax](auto x){return (x/(2*fKmax));});
           std::transform(ky.begin(),ky.end(),trajectories2.begin()+nsamples,[fKmax](auto x){return (x/(2*fKmax));});
            //permute is faster?

            hoNDArray<float> trajectories(2,ADCsamples,Nints_);
            std::vector<size_t> dim_order(trajectories2.get_number_of_dimensions());
            for(int i=0; i<trajectories2.get_number_of_dimensions();i++)
                dim_order[i]=i;

            dim_order[0]=2;
            dim_order[1]=0;
             dim_order[2]=1;
            Gadgetron::permute(trajectories2,trajectories,dim_order);


           
            // GDEBUG("writing trajectories to text  files: /tmp/gadgetron/Trajectories.txt \n");
            // std::ofstream myfile2;
            // myfile2.open ("/tmp/gadgetron/Trajectories.txt",std::ofstream::trunc);
            // myfile2 <<"kx [1/m]          "<<" ky[1/m]         "<<"DCF        "<<'\n';
            // for(auto i=0;i<nsamples*2;i+=2){
            //     myfile2 <<trajectories[i] <<"  "<<trajectories[i+1] <<"  "<<dcf[i/2]<<'\n';
            //     if((i+1)%ADCsamples==0) myfile2 <<"\n\n";
            // }
            // myfile2<<std::flush;
            // myfile2.close();
             
           

            //trajectories.reshape({ADCsamples,Nints_});
            hoNDArray<float> weights(nsamples);
            std::transform(dcf.begin(),dcf.begin()+nsamples,weights.begin(),[](auto x){return x;});
            weights.reshape({ADCsamples,Nints_});

            // std::vector<float> gx,gy;
            // gx=m_SpiralTraj.getGradX();
            // gy=m_SpiralTraj.getGradX();

            // std::ofstream myfile3;
            // myfile3.open ("/tmp/gadgetron/Grad.txt",std::ofstream::trunc);
            // myfile3 <<"gx[mT/m]          "<<" gy[mT/m]         "<<'\n';
            // for(auto i=0;i<gx.size();++i){
            //     myfile3 <<gx[i] <<"  "<<gy[i]<<"  "<<'\n';
            // }
            // myfile3<<std::flush;
            // myfile3.close();


            // if (this->girf_kernel){
            //     base_gradients = correct_gradients(base_gradients,Tsamp_ns_*1e-3,this->girf_sampling_time_us,acq_header.read_dir,acq_header.phase_dir,acq_header.slice_dir);
            //     //Weights should be calculated without GIRF corrections according to Hoge et al 2005
            //     trajectories = calculate_trajectories(base_gradients,sample_time,krmax_);
            // }


            return std::make_pair(std::move(trajectories), std::move(weights));

        }

        
    }
}
