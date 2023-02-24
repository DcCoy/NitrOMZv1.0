function mprint_fig(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default arguments:
A.fig       = - 1;
A.sty       = 'nor';
A.name      = 'figure';
A.for       = 'jpeg';
A.silent    = 0;
% Parse required variables, substituting defaults where necessary
Param = parse_pv_pairs(A, varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Param.for,'jpg')
   Param.for = 'jpeg';
end

if strcmp(Param.name,'figure')&Param.fig~=-1
    Param.name = ['figure' num2str(Param.fig)];
end

if Param.fig~=-1
   figure(Param.fig);
end

if strcmp(Param.sty,'nor') | strcmp(Param.sty,'nor1')  
   set(gcf,'PaperPosition',[0.5 0.5 8.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor2') 
   set(gcf,'PaperPosition',[0.5 0.5 6.0 8.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor2b') 
   set(gcf,'PaperPosition',[0.5 0.5 10.0 14.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor2c') 
   set(gcf,'PaperPosition',[0.5 0.5 6.0 14.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor2d') 
   set(gcf,'PaperPosition',[0.5 0.5 2*6.0 4.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor3') 
   set(gcf,'PaperPosition',[0.5 0.5 6.0 5.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor3b') 
   set(gcf,'PaperPosition',[0.5 0.5 6.0 4.5],'Renderer','Painters');
elseif strcmp(Param.sty,'nor4') 
   set(gcf,'PaperPosition',[0.5 0.5 5.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor4b') 
   set(gcf,'PaperPosition',[0.5 0.5 6.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor5') 
   set(gcf,'PaperPosition',[0.5 0.5 4.0 3.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor6') 
   set(gcf,'PaperPosition',[0.5 0.5 3.0 4.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7') 
   set(gcf,'PaperPosition',[0.5 0.5 10.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7b') 
   set(gcf,'PaperPosition',[0.5 0.5 9.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7c') 
   set(gcf,'PaperPosition',[0.5 0.5 20.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7c1') 
   set(gcf,'PaperPosition',[0.5 0.5 16.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7c2') 
   set(gcf,'PaperPosition',[0.5 0.5 16.0 9.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7d') 
   set(gcf,'PaperPosition',[0.5 0.5 10.0 6.5],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7e') 
   set(gcf,'PaperPosition',[0.5 0.5 12.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7f') 
   set(gcf,'PaperPosition',[0.5 0.5 12.0 9.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor7g') 
   set(gcf,'PaperPosition',[0.5 0.5 14.0 9.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor8') 
   set(gcf,'PaperPosition',[0.5 0.5 9.0 10.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor8b2') 
   set(gcf,'PaperPosition',[0.5 0.5 7.0 7.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor8b') 
   set(gcf,'PaperPosition',[0.5 0.5 9.0 9.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor8c') 
   set(gcf,'PaperPosition',[0.5 0.5 8.0 14.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor9') 
   set(gcf,'PaperPosition',[0.5 0.5 6.0 9.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor10') 
   set(gcf,'PaperPosition',[0.5 0.5 15.0 9.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor11') 
   set(gcf,'PaperPosition',[0.5 0.5 9.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor12') 
   set(gcf,'PaperPosition',[0.5 0.5 11.0 9.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor13') 
   set(gcf,'PaperPosition',[0.5 0.5 9.0 12.0],'Renderer','Painters');
elseif strcmp(Param.sty,'nor13b') 
   set(gcf,'PaperPosition',[0.5 0.5 12 13],'Renderer','Painters');
elseif strcmp(Param.sty,'nor14') 
   set(gcf,'PaperPosition',[0.5 0.5 20.0 12.0],'Renderer','Painters');
elseif strcmp(Param.sty,'sq1') 
   set(gcf,'PaperPosition',[0.5 0.5 6.0 4.8],'Renderer','Painters');
elseif strcmp(Param.sty,'sq2') 
   set(gcf,'PaperPosition',[0.5 0.5 8.0 6.4],'Renderer','Painters');
elseif strcmp(Param.sty,'sq3') 
   set(gcf,'PaperPosition',[0.5 0.5 12.0 10.5],'Renderer','Painters');
elseif strcmp(Param.sty,'sq4') 
   set(gcf,'PaperPosition',[0.5 0.5 10.5 10.5],'Renderer','Painters');
elseif strcmp(Param.sty,'pdf_land') 
   set(gcf,'PaperPosition',[0.2 3.0 8.0 6.0],'Renderer','Painters');
elseif strcmp(Param.sty,'pdf_prof') 
   set(gcf,'PaperPosition',[0.0 0.0 9.0 11.5],'Renderer','Painters');
else 
  error(['invalid format']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if A.silent~=0
   disp(['print -d' Param.for ' ' Param.name ]);
end
eval(['print -d' Param.for ' ' Param.name ]);

