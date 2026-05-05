AbPMg4NPP.m

Main MATLAB function for estimating depth-integrated net primary production (NPP) using the absorption-based primary productivity model (AbPM) framework.

Inputs:
- PAR0   : surface PAR above sea surface
- aph443 : phytoplankton absorption coefficient at 443 nm
- a490   : total absorption coefficient at 490 nm
- bb490  : total backscattering coefficient at 490 nm
- fm     : machine-learning-derived maximum quantum yield of phytoplankton photosynthesis

Output:
- PP     : depth-integrated NPP (mg C m⁻² d⁻¹)
