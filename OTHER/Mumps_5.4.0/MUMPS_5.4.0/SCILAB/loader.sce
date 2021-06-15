path= get_absolute_file_path('loader.sce');
exec(path+"/loader_inc.sce");

functions1 = ["dmumpsc"];
functions2 = ["zmumpsc"];
entrypoint1 = "scidmumps";
entrypoint2 = "scizmumps"; 

addinter(objects,entrypoint1,functions1)
num_interface = floor(funptr("dmumpsc")/100);
intppty(num_interface)

addinter(objects,entrypoint2,functions2)
 num_interface = floor(funptr("zmumpsc")/100);
 intppty(num_interface)

[units,typs,nams]=file();
clear units typs
for k=size(nams,'*'):-1:1 
  l=strindex(nams(k),'loader.sce');
  if l<>[] then
    DIR_SCIMUMPS = part(nams(k),1:l($)-1);
    break
  end
end

DIR_SCIMUMPS_DEM=DIR_SCIMUMPS+ "examples/";

getf(DIR_SCIMUMPS+"initmumps.sci")
getf(DIR_SCIMUMPS+"dmumps.sci")
getf(DIR_SCIMUMPS+"zmumps.sci")

add_help_chapter("Interface to the MUMPS package",path+"Help");


