% %% FILES
ccc;
% fn = 'square2d'; % file root name
% fid.cbt = fopen(sprintf('%s.jou',fn),'w+'); % cubit journal file
% fid.mat = fopen(sprintf('%s.mate',fn),'w+'); % speed mate file
% lab = 0;
% fnc = 0;
% %% DEFINE PROPERTIES
% % _materials_
% 
% % number of materials
% mat.lab = [];
% mat.nm  = 1;
% % material labels
% for i_ = 1:mat.nm
%     mat.lab = [mat.lab;lab+1]; lab = (mat.lab(end)>lab)*mat.lab(end);
% end
% 
% % spectral order
% mat.spo = [4;4];
% % unit weight [kg/m3]
% mat.rho = [1800.;1800.];
% % shear wave velocity [m/s]
% mat.vs = [300;300];
% % lamé coefficients - mu [Pa]
% mat.lmm = mat.vs.^2.*mat.rho;
% % lamé coefficients - nu [1]
% mat.lmn = [0.25;0.25];
% % P-wave modulus [Pa]
% mat.lmp = 2*mat.lmm.*(1-mat.lmn)./(1-2*mat.lmn);
% % pressure wave velocity [m/s]
% mat.vp  = sqrt(mat.lmp./mat.rho);
% % yield strain
% mat.eyl = [1e-6;1e-6];
% % yield strength
% mat.syl = sqrt(3).*mat.lmm.*mat.eyl;
% % tau_lim/tau_yld ratio
% mat.ttr = [10;10];
% % tau limit
% mat.tlm = mat.ttr.*mat.syl;
% % ultimate strength (clays)
% mat.su = mat.tlm./sqrt(3);
% % hardening parameters
% mat.ckn = 2.*mat.lmm.*(1+mat.lmn);
% mat.kkn = mat.ckn./(mat.tlm-mat.syl);
% % elastic shear modulus [Pa]
% mat.lme = mat.ckn;
% % lamé coefficients - lambda [Pa]
% mat.lml = mat.lme.*mat.lmn./(1+mat.lmn)./(1-2*mat.lmn);
% % isotropic hardening
% mat.bis = [0;0];
% mat.ris = [0;0];
% % nonlinear material
% mat.nln = 1;
% % viscous gamma
% mat.vgm = [0.;0.];
% 
% %%
% % _geometry_
% 
% % number of surfaces
% geo.nq  = mat.nm;
% % surface labels
% geo.lab.srf = mat.lab+0:6:6*geo.nq-1;
% % surface side lengths
% geo.xdm = [1.,1.];
% geo.ydm = [1.;1.];
% 
% %%
% % _mesh size_
% 
% % number of elements per material
% msh.ne  = [1;1];
% % mesh sizes
% msh.xms = [1;1];
% msh.yms = [1;1];
% 
% %%
% % _absorbing bc_
% 
% % number of abc conditions
% % abc.nc = 0;
% % % abc labels
% % abc(1).lab = lab+1; lab = (abc(1).lab(end)>lab)*abc(1).lab(end);
% % % abc curves
% % abc(1).crv{1} = [1;3];
% % for i_ = 1:abc.nc-1
% %     abc(1).crv{1} = [abc(1).crv{1};abc(1).crv{1}+12*i_];
% % end
% 
% %%
% % _dirichlet bc_
% % number of drx conditions
% drx.nc = 3;
% % drx labels
% drx(1).lab = lab+1; lab = (drx(1).lab(end)>lab)*drx(1).lab(end);
% % drx curves
% drx(1).crv{1} = 1;
% % drx functions
% drx(1).fnc = fnc+1; fnc = (drx(1).fnc(end)>fnc)*drx(1).fnc(end);
% % drx type
% drx(1).typ = 12;
% % drx values
% drx(1).vls = [0;1];
% % drx parameters
% drx(1).prm =  [10,10,1];
% 
% % drx labels
% drx(2).lab = lab+1; lab = (drx(2).lab(end)>lab)*drx(2).lab(end);
% % drx curves
% drx(2).crv{1} = 3;
% % drx functions
% drx(2).fnc = fnc+1; fnc = (drx(2).fnc(end)>fnc)*drx(2).fnc(end);
% % drx type
% drx(2).typ = 12;
% % drx values
% drx(2).vls = [1;0];
% % drx parameters
% drx(2).prm =  [10,10,1];
% 
% 
% % drx labels
% drx(3).lab = lab+1; lab = (drx(3).lab(end)>lab)*drx(3).lab(end);
% % drx curves
% drx(3).crv{1} = 4;
% % drx functions
% drx(3).fnc = fnc+1; fnc = (drx(3).fnc(end)>fnc)*drx(3).fnc(end);
% % drx type
% drx(3).typ = 0;
% % drx values
% drx(3).vls = [0;0];
% % drx parameters
% drx(3).prm = 0;
% 
% % number of dry conditions
% dry.nc = 2;
% % dry labels
% dry(1).lab = drx(3).lab;%lab+1; lab = (dry(1).lab(end)>lab)*dry(1).lab(end);
% % dry curves
% dry(1).crv{1} = 4;
% % dry functions
% dry(1).fnc = fnc+1; fnc = (dry(1).fnc(end)>fnc)*dry(1).fnc(end);
% % dry type
% dry(1).typ = 0;
% % dry values
% dry(1).vls = [0;0];
% % dry parameters
% dry(1).prm = 0;
% 
% % dry labels
% dry(2).lab = lab+1; lab = (dry(2).lab(end)>lab)*dry(2).lab(end);
% % dry curves
% dry(2).crv{1} = 2;
% % dry functions
% dry(2).fnc = fnc+1; fnc = (dry(2).fnc(end)>fnc)*dry(2).fnc(end);
% % dry type
% dry(2).typ = 0;
% % dry values
% dry(2).vls = [0;0];
% % dry parameters
% dry(2).prm = 0;
% 
% %%
% % _neumann bc_
% % number of nex conditions
% % nex.nc = 0;
% % % nex labels
% % nex(1).lab = lab+1; lab = (nex(1).lab(end)>lab)*nex(1).lab(end);
% % % nex curves
% % nex(1).crv(1) = 2+12*(mat.nm-1);
% % % nex functions
% % nex(1).fnc = fnc+1; fnc = (nex(1).fnc(end)>fnc)*nex(1).fnc(end);
% % % nex type
% % nex(1).typ = 12;
% % % nex values
% % nex(1).vls = [1;1];
% % % nex parameters
% % nex(1).prm = [1000,10,1];
% 
% %%
% % _monitor_
% 
% % % number of monitors
% % mon.nm = 2;
% % mon(1).lab = lab+1; lab = (mon(1).lab(end)>lab)*mon(1).lab(end);
% % mon(1).x = [geo.xdm(1)/2;geo.xdm(1)/2];
% % mon(1).y = [0;geo.ydm(2)];
% 
% %% CUBIT JOURNAL FILE
% % _create main brick_
% 
% fprintf(fid.cbt,'reset\n');
% xof = 0;
% yof = 0;
% for i_ = 1:mat.nm
%     fprintf(fid.cbt,'brick x %f y %f\n',geo.xdm(i_),geo.ydm(i_));
%     fprintf(fid.cbt,'move surface %u x %f y %f include_merged\n',geo.lab.srf(i_),...
%         geo.xdm(i_)/2+xof,geo.ydm(i_)/2+yof);
%     yof = yof+geo.ydm(i_);
% end
% 
% %%
% % _mesh curves/surfaces_
% mc = 0;
% for j_ = 1:mat.nm
%     for i_ = 1:2
%         fprintf(fid.cbt,'curve %u interval %.1f\n',2*i_-1+mc,...
%             msh.xms(i_));
%         fprintf(fid.cbt,'curve %u interval %.1f\n',2*i_+mc,...
%             msh.yms(i_));
%         fprintf(fid.cbt,'mesh curve %u %u\n',2*i_-1+mc,...
%             2*i_+mc);
%     end
%     mc = mc+12;
%     fprintf(fid.cbt,'mesh surface %u\n',geo.lab.srf(j_));
%     fprintf(fid.cbt,'block %u surface %u\n',mat.lab(j_),geo.lab.srf(j_));
%     fprintf(fid.cbt,'block %u element type quad\n',mat.lab(j_));
% end
% fprintf(fid.cbt,'merge all\n');
% 
% %%
% % _boundary conditions_
% % fprintf(fid.cbt,['block %u curve',repmat(' %u',1,numel(abc(1).crv{:})),'\n'],...
% %     abc(1).lab,abc(1).crv{:});
% fprintf(fid.cbt,'block %u curve %u \n',drx(1).lab,drx(1).crv{:});
% fprintf(fid.cbt,'block %u curve %u \n',drx(2).lab,drx(2).crv{:});
% fprintf(fid.cbt,'block %u curve %u \n',drx(3).lab,drx(3).crv{:});
% % fprintf(fid.cbt,'block %u curve %u \n',dry(1).lab,dry(1).crv{:});
% fprintf(fid.cbt,'block %u curve %u \n',dry(2).lab,dry(2).crv{:});
% % fprintf(fid.cbt,'block %u curve %u \n',nex.lab,nex.crv);
% 
% %% LS INPUT FILE
% %
% % fprintf(fid.cbt,'create vertex %f %f 0.5\n',mon.x(1),mon.y(1));
% % fprintf(fid.cbt,'create vertex %f %f 0.5\n',mon.x(2),mon.y(2));
% % fprintf(fid.cbt,'create curve vertex 9 to 10\n');
% % fprintf(fid.cbt,'curve 13 interval 1\n');
% % fprintf(fid.cbt,'mesh curve 13\n');
% % fprintf(fid.cbt,'block %u curve 13\n',mon.lab);
% 
% %%
% % _export mesh_
% 
% fprintf(fid.cbt,'set large exodus file off\n');
% fprintf(fid.cbt,['export mesh "%s.e" dimension 2 block',...
%     repmat(' %u',1,lab),' overwrite\n'],fn,1:lab);
% % fprintf(fid.cbt,'undo group begin\n');
% % fprintf(fid.cbt,'set large exodus file off\n');
% % fprintf(fid.cbt,'export mesh "LS%s.e" dimension 2 block %u overwrite\n',fn,mon.lab);
% % fprintf(fid.cbt,'undo group end\n');
% % fprintf(fid.cbt,'Exit');
% fclose(fid.cbt);


%% CREATE AND CONVERT MESH
str = ['! wine /media/filippo/OS/Program\ Files\ \(x86\)/CUBIT\ ',...
    '13.0/bin/cubit.exe -nographics -input /media/filippo/Data/Filippo/',...
    'PHD_passing_through/MATLAB/m_polimi/square2d.jou'];
eval(str);

str = ['! ncdump /media/filippo/Data/Filippo/PHD_heavyweight/SPEED2D/INPUTS/',...
    'plastic/square2d.e > /media/filippo/Data/Filippo/PHD_heavyweight/SPEED2D/INPUTS/',...
    'plastic/square2d.txt'];
eval(str);
exo2mesh_2D(['/media/filippo/Data/Filippo/PHD_heavyweight/SPEED2D/INPUTS/',...
    'plastic/'],'square2d');
% str = sprintf('! ncdump LS%s.e > LS%s.txt',fn,fn);
% eval(str);
% exo2mesh_3Dtria('./',sprintf('LS%s',fn));

% %% MATE FILE
% % _material properties_
% if mat.nln
%     for i_ = 1:mat.nm
%         fprintf(fid.mat,'MATE%6u%6u%10.1f%12.4e%12.4e%12.1e%12.4e%12.4e%12.4e%6.1f%6.1f\n\n',...
%             mat.lab(i_),mat.spo(i_),mat.rho(i_),mat.lml(i_),mat.lmm(i_),...
%             mat.vgm(i_),mat.syl(i_),mat.ckn(i_),mat.kkn(i_),mat.bis(i_),...
%             mat.ris(i_));
%     end
% else
%     for i_ = 1:mat.nm
%         fprintf(fid.mat,'MATE%10u%10u%15.1f%15.4e%15.4e%15.4e\n\n',...
%             mat.lab(i_),mat.spo(i_),mat.rho(i_),mat.lml(i_),mat.lmm(i_),...
%             mat.vgm(i_));
%     end
% end
% %%
% % _absorbing bc_
% % fprintf(fid.mat,'ABSO %u\n\n',abc(1).lab);
% %%
% % _X dirichlet bc_
% fprintf(fid.mat,'DIRX %u %u %f %f\n',drx(1).lab,drx(1).fnc,...
%     drx(1).vls(1),drx(1).vls(2));
% fprintf(fid.mat,'FUNC %u %u %f %f %f\n\n',drx(1).fnc,drx(1).typ,drx(1).prm);
% %
% fprintf(fid.mat,'DIRX %u %u %f %f\n',drx(2).lab,drx(2).fnc,...
%     drx(2).vls(1),drx(2).vls(2));
% fprintf(fid.mat,'FUNC %u %u %f %f %f\n\n',drx(2).fnc,drx(2).typ,drx(2).prm);
% %
% fprintf(fid.mat,'DIRX %u %u %f %f\n',drx(3).lab,drx(3).fnc,...
%     drx(3).vls(1),drx(3).vls(2));
% fprintf(fid.mat,'FUNC %u %u %f\n\n',drx(3).fnc,drx(3).typ,drx(3).prm);
% %%
% % _Y dirichlet bc_
% fprintf(fid.mat,'DIRY %u %u %f %f\n',dry(1).lab,dry(1).fnc,...
%     dry(1).vls(1),dry(1).vls(2));
% fprintf(fid.mat,'FUNC %u %u %f\n\n',dry(1).fnc,dry(1).typ,dry(1).prm);
% %
% fprintf(fid.mat,'DIRY %u %u %f %f\n',dry(2).lab,dry(2).fnc,...
%     dry(2).vls(1),dry(2).vls(2));
% fprintf(fid.mat,'FUNC %u %u %f\n\n',dry(2).fnc,dry(2).typ,dry(2).prm);
% % %%
% % % _neumann bc_
% % fprintf(fid.mat,'NEUX %u %u %f %f\n',nex.lab,nex.fnc,...
% %     nex.vls(1),nex.vls(2));
% % fprintf(fid.mat,'FUNC %u %u %f %f %f\n\n',nex.fnc,nex.typ,...
% %     nex.prm(1),nex.prm(2),nex.prm(3));
% % fclose(fid.mat);