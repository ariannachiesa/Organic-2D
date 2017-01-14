function [omesh] = mesh2d_reduction(device)

    omesh = device.msh;
    omesh.p = device.msh.p(:,device.scnodes);
    omesh.t = device.msh.t(:,! device.insulator);
    omesh.e = device.msh.e(:,device.msh.e(end,:)==1);
    omesh.wjacdet = device.msh.wjacdet(:,! device.insulator);
    omesh.area = device.msh.area(! device.insulator,:);
    omesh.shg = device.msh.shg(:,:,! device.insulator);
    
endfunction