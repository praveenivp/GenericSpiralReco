# Generic Gadgetron reconstruction chain for peSpiral Sequence


## Dependentcies
    - Gadgetron
    - BART

## Building
```
cmake -DCMAKE_BUILD_TYPE=Debug -GNinja ..
ninja
```

## bart troubleshooting
make sure you enable `calc_csm` and `do_pics` flag
go the output folder and execute the following
```
bart pics -i 10  -p DCF -t Traj coil_data sensCFL pics_reco

```

## getting ISMRMRD files for offline reco
```
siemens_to_ismrmrd -x parameter_maps/IsmrmrdParameterMap_peSpiral_VE11E.xsl -m parameter_maps/IsmrmrdParameterMap_peSpiral_VE11E.xml -f peSpiral.dat -z 1 -o data_z1.h5
ismrmrdviewer data_z1.h5 
gadgetron_ismrmrd_client -f testDaat/data_z1.h5 -C config/piv_GadgetronSpiral.xml 
```

# Tunneling to your own Gadgetron
Turn off tunneling in `gadgetron.ini` file and make your own tunnel. 
```
ssh -L 192.168.2.1:9010:localhost:9002 gadgetron@10.41.116.43
```

# deploying in docker-conda environment for online reconstruction
mount the folder and start the docker container
```
scp -r <folder_path> username@192.168.0.1:~
docker run --rm -v <GenericSpiralReco_PATH>:/code <container_image_id>
docker exec -ti <container_id> bash
```
Inside the container shell, make sure the linker (ld) points to the `ld` inside conda by `which ld`. Otherwise, pthread library path was incorrect(`ldconfig -p |grep pthread`)
```
conda activate gadgetron
export PATH="${CONDA_PREFIX}/x86_64-conda-linux-gnu/bin":$PATH
 ```
install missing packages and compile. If sucessful, commit the `lib.so` to [gadgetron-gagets]() repo and restart gadgetron container listening to port 9002 vai yacht.
```
apt update
apt install cmake libarmadillo-dev
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} ..
make -j8
scp <libfile.so> gadgetron@10.41.116.43:/home/gadgetron/Documents/packages/gadgetron-gadgets/lib
```



## Todo  
- [x]  Accumulate data
- [x]  Add Trajectory information to Ismrmrd Reconbit
- [x]  Add nice NODE_PROPERTY for tunning gadgets with xml config.
- [x]  import and export functions for CFL files and h0NDArrays.
- [x]  Calculte coil maps from reference data using BART
- [x]  gridding reoconstruction
- [x]  coil combination with calculated sensitivity
- [x]  CG-SENSE reconstruction
- []  remove CFL_IO2.xx 
- []  multi-slice reco fix
- []  GIRF correction
- []  coil maps calculation not starting immediately
- []  Mosiacing doesn't work with reference image

