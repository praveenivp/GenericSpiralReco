# developing and testing in a remote apptainer container
The conda environment proposed for development is not optimal. let's see whether we can develop and test new C++ gadgets in a remote apptainer container.

## building container
first build the apptiner container 
`appatainer build gadgetron.sif gadgetron.recipe`

you might need `ismrmrdviewer` from `mrdviewer.recipe` for seeing the images.

## Local container

need to test

## remote container

### setting up remote explorer in vscode
add the following lines to your ssh config file
```
Host RemoteGagetron
	RemoteCommand apptainer shell -B /ptmp /ptmp/pvalsala/MyContainers/gadgetron.sif
	RequestTTY yes
	User pvalsala
    HostName nyx.hpc.kyb.local
	ForwardX11 yes
	ForwardX11Trusted yes
```
`RemoteGagetron` should appear in your remote explorer

## building the library
Now you can open this repo cloned on the remote. Opening terminal will give you a terminal inside the container.
```
mkcd build
cmake .. 
make
```

## testing

### simple testing
	open another terminal in vscode and launch gadgetron.

### testing in GPU
```
salloc-a-node --time 02:00:00 --partition gpu
```
keep that open another terminal and ssh into it with portforwarding and start gadgetron
```
ssh -Y -L 9004:localhost:9002 gpu-x
apptainer exec -B /ptmp  --nv /ptmp/pvalsala/MyContainers/gadgetron.sif bash
cd /so/file/location
gadgetron
```
now you can send data to gadgetron from login node to port 9004

```
gadgetron_ismrmrd_client -f testData/raw3D_R3.h5 -C config/piv_GadgetronSpiral_fast.xml -o im_2.h5 -p 9004
```


## troubleshooting
- `Ctrl+Shift+P` : search for `kill vs code server` in vscode
- delete `/home/user/.vscode-server` folder
- `ps aux | grep -E ${USER}.*.vscode-server | awk '{print $2}' |xargs kill -9` to kill all vscode-servers running with your user.


## references
-https://code.visualstudio.com/docs/devcontainers/containers