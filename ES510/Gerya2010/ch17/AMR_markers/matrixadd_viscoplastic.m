% Function matrixad_test4()
% This function compose matrix and right parts
% for Poisson and then for continuity and Stokes equations
% Function return solutions for Poisson SP() and Stokes S()
function[S, exx, eyy, exy, sxx, syy, sxy, pr]=matrixadd_viscoplastic(varnum,celnum,nodnum,celdot,celnod,celvar,celres,celeta,celrho,nodx,nody,nodcel,nodeta,nodrho,bupper,blower,bleft,bright,bpres,g,pscale,gx,gy)

% Define global matrices for Stokes
L=sparse(varnum,varnum);
R=zeros(varnum,1);
% Processed Y/N index
varyn=zeros(varnum,1);

% Composing equations
for ci=1:1:celnum
    if(celdot(ci,1)==0)
        %
        %       |ciA1 |ciA2 | 
        %   ----+----vy2----+----
        %  ciL1 |  ci       | ciR1
        %   ---vx1   P5    vx4---
        %  ciL2 |           | ciR2
        %   ----+----vy3----+----
        %       |ciB1 |ciB2 |
        % Indexes of cells
        % Left
        ciL1=nodcel(celnod(ci,1),2); 
        ciL2=nodcel(celnod(ci,2),1);
        % Right
        ciR1=nodcel(celnod(ci,3),4);
        ciR2=nodcel(celnod(ci,4),3);
        % Above
        ciA1=nodcel(celnod(ci,1),3); 
        ciA2=nodcel(celnod(ci,3),1); 
        % Below
        ciB1=nodcel(celnod(ci,2),4);
        ciB2=nodcel(celnod(ci,4),2);
        
        % Equation for P5------------------------
        % Index for P5
        kp=celvar(ci,5);
        % Boundary condition
        if(nodx(celnod(ci,1))==g.xmin && nody(celnod(ci,1))==g.ymin)
                % Upper-left cell
                L(kp,kp)=1*pscale;
                R(kp)= bpres;
       else
            % Continuity equation
            %     vy2
            % vx1  P5   vx4
            %     vy3
            kfx=1/(nodx(celnod(ci,3))-nodx(celnod(ci,1)));
            kfy=1/(nody(celnod(ci,2))-nody(celnod(ci,1)));
            L(kp,celvar(ci,1))=-kfx; % vx1
            L(kp,celvar(ci,2))=-kfy; % vy2
            L(kp,celvar(ci,3))=+kfy; % vy3
            L(kp,celvar(ci,4))=+kfx; % vx4
            R(kp)=0;
        end
        % Mark processed equation
        varyn(kp)=1;
        % End Equation for P5------------------------
        
        
        % Equation for vx1------------------------
        % Index for vx1
        kx=celvar(ci,1);
        % Check processed variable y/n
        if(varyn(kx)==0)
            % Boundary condition
            if(nodx(celnod(ci,1))==g.xmin)
                L(kx,kx)=1;
                R(kx)=bleft(1);
                % Mark processed equation
                varyn(kx)=1;
            else
                % Solving x-Stokes
                % dSIGxx/dx-dP/dx + dSIGxy/dy - dt/2*gx*(vx*dRHO/dx+vy*dRHO/dy)= -RHO*gx
                % Cell/Subcells -> Cell
                if(ciL1>0 && celres(ciL1)>=celres(ci))
                    % Compose x-Stokes equation 
                    % dSIGxx/dx-dP/dx + dSIGxy/dy = -RHO*gx
                    %            SIGxy1
                    % SIG'xx1,P1    vx1    SIG'xx2,P2    vx4
                    %            SIGxy2
                    % Compose matrix
                    % dSIGxx/dx-dP/dx
                    % x-distance between SIGxx-nodes
                    dx=abs((nodx(celnod(ci,3))-nodx(celnod(ciL1,1)))/2); 
                    % Add SIGxx2=SIG'xx2-P2
                    [L]=sigmaxx_test(kx,ci,nodx,celnod,celvar,celeta,1/dx,pscale/dx,L);
                    % Add SIGxx1=SIG'xx1-P1
                    if(celres(ciL1)==celres(ci))
                        % No change in resolution
                        [L]=sigmaxx_test(kx,ciL1,nodx,celnod,celvar,celeta,-1/dx,-pscale/dx,L);
                    else
                        % Increase in resolution
                        % Upper left cell
                        [L]=sigmaxx_test(kx,ciL1,nodx,celnod,celvar,celeta,-0.5/dx,-0.5*pscale/dx,L);                        % P1
                        % Lower left cell
                        [L]=sigmaxx_test(kx,ciL2,nodx,celnod,celvar,celeta,-0.5/dx,-0.5*pscale/dx,L);                        % P1
                    end
                    % dSIGxy/dy
                    % y-distance between SIGxy-nodes
                    dy=abs(nody(celnod(ci,2))-nody(celnod(ci,1))); 
                    % Add SIGxy1
                    if(ciA1==0 && ciA2>0)
                        % vx between two subcells
                        %   |             |
                        % -ni1-----+-----ni2-
                        %   |      |      |       
                        %         vx1
                        % Left top node
                        [L,R]=sigmaxy_test(kx,celnod(ciL1,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,R,bupper,blower,bleft,bright,g);
                        % Right top node
                        [L,R]=sigmaxy_test(kx,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vx in standard position
                        %          |
                        %   ------ni1------
                        %          |             
                        %         vx1
                        % Top node
                        [L,R]=sigmaxy_test(kx,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dy,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Add SIGxy2
                    if(ciB1==0 && ciB2>0)
                        % vx between two subcells
                        %         vx1
                        %   |      |      |       
                        % -ni1-----+-----ni2-
                        %   |             |
                        % Left bottom node
                        [L,R]=sigmaxy_test(kx,celnod(ciL2,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,R,bupper,blower,bleft,bright,g);
                        % Right bottom node
                        [L,R]=sigmaxy_test(kx,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vx in standard position
                        %         vx1
                        %          |
                        %   ------ni1------
                        %          |             
                        % Top node
                        [L,R]=sigmaxy_test(kx,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dy,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Right Part -RHO*gx
%                     R(kx)=R(kx)-gx*(nodrho(celnod(ci,1))+nodrho(celnod(ci,2)))/2;
                    R(kx)=R(kx)-gx*celrho(ci,1);
                    % Mark processed equation
                    varyn(kx)=1;
                else
                    % dvx/dy stress gradient balance
                    dy=abs(nody(celnod(ci,2))-nody(celnod(ci,1)))/2; 
                    % Upper subcel 
                    if(ciL1>0)
%                         L(kx,kx)=-(nodeta(celnod(ciL1,3))+nodeta(celnod(ciL1,4)))/2/dy; % vx1
%                         L(kx,celvar(ciL1,4))=(nodeta(celnod(ciL1,3))+nodeta(celnod(ciL1,4)))/2/dy; % vx4
                        L(kx,kx)=-nodeta(celnod(ci,2))/dy; % vx1
                        L(kx,celvar(ciL1,4))=nodeta(celnod(ci,2))/dy; % vx4
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciL1,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciL1,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    % Lower subcel
                    else
%                         L(kx,kx)=(nodeta(celnod(ciL2,3))+nodeta(celnod(ciL2,4)))/2/dy; % vx1
%                         L(kx,celvar(ciL2,4))=-(nodeta(celnod(ciL2,3))+nodeta(celnod(ciL2,4)))/2/dy; % vx4
                        L(kx,kx)=nodeta(celnod(ci,1))/dy; % vx1
                        L(kx,celvar(ciL2,4))=-nodeta(celnod(ci,1))/dy; % vx4
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciL2,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciL2,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Mark processed equation
                    varyn(kx)=1;
                end
            end
        end
        % End Equation for vx1------------------------
        
        
        % Equation for vx4------------------------
        % Index for vx4
        kx=celvar(ci,4);
        % Check processed variable y/n
        if(varyn(kx)==0)
            % Boundary condition
            if(nodx(celnod(ci,3))==g.xmax)
                L(kx,kx)=1;
                R(kx)= bright(1);
                % Mark processed equation
                varyn(kx)=1;
            else
                % Solving x-Stokes
                % dSIGxx/dx-dP/dx + dSIGxy/dy - dt/2*gx*(vx*dRHO/dx+vy*dRHO/dy)= -RHO*gx
                %  Cell->Cell/Subcells    
                if(ciR1>0 && celres(ciR1)>=celres(ci))
                    % Compose x-Stokes equation 
                    % dSIGxx/dx-dP/dx + dSIGxy/dy = -RHO*gx
                    %                    SIGxy1
                    % vx1   SIG'xx1,P1    vx4    SIG'xx2,P2    
                    %                    SIGxy2
                    % Compose matrix
                    % dSIGxx/dx-dP/dx
                    % x-distance between SIGxx-nodes
                    dx=abs((nodx(celnod(ciR1,3))-nodx(celnod(ci,1)))/2); 
                    % Add SIGxx1=SIG'xx1-P1
                    [L]=sigmaxx_test(kx,ci,nodx,celnod,celvar,celeta,-1/dx,-pscale/dx,L);
                    % Add SIGxx2=SIG'xx2-P2
                    if(celres(ciR1)==celres(ci))
                        % No change in resolution
                        [L]=sigmaxx_test(kx,ciR1,nodx,celnod,celvar,celeta,1/dx,pscale/dx,L);
                    else
                        % Increase in resolution
                        % Upper right cell
                        [L]=sigmaxx_test(kx,ciR1,nodx,celnod,celvar,celeta,0.5/dx,0.5*pscale/dx,L);                        % P1
                        % Lower right cell
                        [L]=sigmaxx_test(kx,ciR2,nodx,celnod,celvar,celeta,0.5/dx,0.5*pscale/dx,L);                        % P1
                    end
                    % dSIGxy/dy
                    % y-distance between SIGxy-nodes
                    dy=abs(nody(celnod(ci,4))-nody(celnod(ci,3))); 
                    % Add SIGxy1
                    if(ciA1>0 && ciA2==0)
                        % vx between two subcells
                        %   |     ciA1    |
                        % -ni1-----+-----ni2-
                        %   |  ci  |      |       
                        %         vx4
                        % Left top node
                        [L,R]=sigmaxy_test(kx,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,R,bupper,blower,bleft,bright,g);
                        % Right top node
                        [L,R]=sigmaxy_test(kx,celnod(ciR1,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dy,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vx in standard position
                        %          |
                        %   ------ni1------
                        %          |             
                        %         vx4
                        % Top node
                        [L,R]=sigmaxy_test(kx,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dy,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Add SIGxy2
                    if(ciB1>0 && ciB2==0)
                        % vx between two subcells
                        %         vx4
                        %   |  ci  |      |       
                        % -ni1-----+-----ni2-
                        %   |     ciB1    |
                        % Left bottom node
                        [L,R]=sigmaxy_test(kx,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,R,bupper,blower,bleft,bright,g);
                        % Right bottom node
                        [L,R]=sigmaxy_test(kx,celnod(ciR2,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dy,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vx in standard position
                        %         vx1
                        %          |
                        %   ------ni1------
                        %          |             
                        % Bottom node
                        [L,R]=sigmaxy_test(kx,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dy,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Right Part
%                     R(kx)=R(kx)-gx*(nodrho(celnod(ci,3))+nodrho(celnod(ci,4)))/2;
                    R(kx)=R(kx)-gx*celrho(ci,4);
                    % Mark processed equation
                    varyn(kx)=1;
                else
                    % dvx/dy stress gradient balance
                    dy=abs(nody(celnod(ci,2))-nody(celnod(ci,1)))/2; 
                    % Upper subcel 
                    if(ciR1>0)
%                         L(kx,kx)=-(nodeta(celnod(ciR1,1))+nodeta(celnod(ciR1,2)))/2/dy; % vx4
%                         L(kx,celvar(ciR1,1))=(nodeta(celnod(ciR1,1))+nodeta(celnod(ciR1,2)))/2/dy; % vx1
                        L(kx,kx)=-nodeta(celnod(ci,4))/dy; % vx4
                        L(kx,celvar(ciR1,1))=nodeta(celnod(ci,4))/dy; % vx1
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciR1,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciR1,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    % Lower subcel
                    else
%                         L(kx,kx)=(nodeta(celnod(ciR2,1))+nodeta(celnod(ciR2,2)))/2/dy; % vx4
%                         L(kx,celvar(ciR2,1))=-(nodeta(celnod(ciR2,1))+nodeta(celnod(ciR2,2)))/2/dy; % vx1
                        L(kx,kx)=nodeta(celnod(ci,3))/dy; % vx4
                        L(kx,celvar(ciR2,1))=-nodeta(celnod(ci,3))/dy; % vx1
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciR2,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vx_test(kx,celnod(ciR2,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Mark processed equation
                    varyn(kx)=1;
                end
            end
        end
        % End Equation for vx4------------------------
        
        
        % Equation for vy2------------------------
        % Index for vy2
        ky=celvar(ci,2);
        % Check processed variable y/n
        if(varyn(ky)==0)
            % Boundary condition
            if(nody(celnod(ci,1))==g.ymin)
                L(ky,ky)=1;
                R(ky)= bupper(1);
                % Mark processed equation
                varyn(ky)=1;
            else
                % Solving y-Stokes
                % dSIGyy/dy-dP/dy + dSIGxy/dx - dt/2*gy*(vx*dRHO/dx+vy*dRHO/dy) = -RHO*gy
                % Cell/Subcells -> Cell   
                if(ciA1>0 && celres(ciA1)>=celres(ci))
                    % Compose y-Stokes equation 
                    % dSIGyy/dy-dP/dy + dSIGxy/dx = -RHO*gy
                    %           SIGyy1,P1
                    %   SIGxy1     vy2    SIG'xy2
                    %           SIGyy2,P2
                    %              vy3
                    % Compose matrix
                    % dSIGyy/dy-dP/dy
                    % y-distance between SIGyy-nodes
                    dy=abs((nody(celnod(ci,2))-nody(celnod(ciA1,1)))/2); 
                    % Add SIGyy2=SIG'yy2-P2
                    [L]=sigmayy_test(ky,ci,nody,celnod,celvar,celeta,1/dy,pscale/dy,L);
                    % Add SIGyy1=SIG'yy1-P1
                    if(celres(ciA1)==celres(ci))
                        % No change in resolution
                        [L]=sigmayy_test(ky,ciA1,nody,celnod,celvar,celeta,-1/dy,-pscale/dy,L);
                    else
                        % Increase in resolution
                        % Upper left cell
                        [L]=sigmayy_test(ky,ciA1,nody,celnod,celvar,celeta,-0.5/dy,-0.5*pscale/dy,L);                        % P1
                        % Upper right cell
                        [L]=sigmayy_test(ky,ciA2,nody,celnod,celvar,celeta,-0.5/dy,-0.5*pscale/dy,L);                        % P1
                   end
                    % dSIGxy/dx
                    % x-distance between SIGxy-nodes
                    dx=abs(nodx(celnod(ci,3))-nodx(celnod(ci,1))); 
                   % Add SIGxy1
                    if(ciL1==0 && ciL2>0)
                        % vy between two subcells
                        %       |
                        %------ni1------
                        %       | ciA1
                        %  ciL2 +--vy2
                        %       | ci
                        %------ni2------
                        %       |
                        % Upper left node
                        [L,R]=sigmaxy_test(ky,celnod(ciA1,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,R,bupper,blower,bleft,bright,g);
                        % Lower left node
                        [L,R]=sigmaxy_test(ky,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vy in standard position
                        %          |
                        %   ------ni1---vy2
                        %          |  ci           
                        % Left node
                        [L,R]=sigmaxy_test(ky,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dx,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Add SIGxy2
                    if(ciR1==0 && ciR2>0)
                        % vy between two subcells
                        %         |
                        %  ------ni1------
                        %    ciA2 | 
                        %  --vy2--+ ciR2
                        %     ci  | 
                        %  ------ni2------
                        %         |
                        % Left bottom node
                        [L,R]=sigmaxy_test(ky,celnod(ciA2,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,R,bupper,blower,bleft,bright,g);
                        % Right bottom node
                        [L,R]=sigmaxy_test(ky,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vy in standard position
                        %          |
                        %   vy2---ni1---
                        %     ci   |            
                        % Top node
                        [L,R]=sigmaxy_test(ky,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dx,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Right Part
%                     R(ky)=R(ky)-gy*(nodrho(celnod(ci,1))+nodrho(celnod(ci,3)))/2;
                    R(ky)=R(ky)-gy*celrho(ci,2);
                    % Mark processed equation
                    varyn(ky)=1;
                else
                    % dvy/dx stress gradient balance
                    dx=abs(nodx(celnod(ci,3))-nodx(celnod(ci,1)))/2; 
                    % Left subcel
                    if(ciA1>0)
%                         L(ky,ky)=-(nodeta(celnod(ciA1,2))+nodeta(celnod(ciA1,4)))/2/dx; % vy2
%                         L(ky,celvar(ciA1,3))=(nodeta(celnod(ciA1,2))+nodeta(celnod(ciA1,4)))/2/dx; % vy3
                        L(ky,ky)=-nodeta(celnod(ci,3))/dx; % vy2
                        L(ky,celvar(ciA1,3))=nodeta(celnod(ci,3))/dx; % vy3
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciA1,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciA1,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    % Right subcel
                    else
%                         L(ky,ky)=(nodeta(celnod(ciA2,2))+nodeta(celnod(ciA2,4)))/2/dx; % vy2
%                         L(ky,celvar(ciA2,3))=-(nodeta(celnod(ciA2,2))+nodeta(celnod(ciA2,4)))/2/dx; % vy3
                        L(ky,ky)=nodeta(celnod(ci,1))/dx; % vy2
                        L(ky,celvar(ciA2,3))=-nodeta(celnod(ci,1))/dx; % vy3
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciA2,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciA2,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Mark processed equation
                    varyn(ky)=1;
                end
            end
        end
        % End Equation for vy2------------------------
                
        
        
        % Equation for vy3------------------------
        % Index for vy3
        ky=celvar(ci,3);
        % Check processed variable y/n
        if(varyn(ky)==0)
            % Boundary condition
            if(nody(celnod(ci,4))==g.ymax)
                L(ky,ky)=1;
                R(ky)= blower(1);
                % Mark processed equation
                varyn(ky)=1;
            else
                % Solving y-Stokes
                % dSIGyy/dy-dP/dy + dSIGxy/dx - dt/2*gy*(vx*dRHO/dx+vy*dRHO/dy) = -RHO*gy
                % Cell/Subcells -> Cell   
                if(ciB1>0 && celres(ciB1)>=celres(ci))
                    % Compose y-Stokes equation 
                    % dSIGyy/dy-dP/dy + dSIGxy/dx = -RHO*gy
                    %              vy2
                    %           SIGyy1,P1
                    %   SIGxy1     vy3    SIG'xy2
                    %           SIGyy2,P2
                    % Compose matrix
                    % dSIGyy/dy-dP/dy
                    % y-distance between SIGyy-nodes
                    dy=abs((nody(celnod(ciB1,2))-nody(celnod(ci,1)))/2); 
                    % Add SIGyy1=SIG'yy1-P1
                    [L]=sigmayy_test(ky,ci,nody,celnod,celvar,celeta,-1/dy,-pscale/dy,L);
                    % Add SIGyy2=SIG'yy2-P2
                    if(celres(ciB1)==celres(ci))
                        % No change in resolution
                        [L]=sigmayy_test(ky,ciB1,nody,celnod,celvar,celeta,1/dy,pscale/dy,L);
                        % Compute gy=-dFI/dy
                   else
                        % Increase in resolution
                        % Upper left cell
                        [L]=sigmayy_test(ky,ciB1,nody,celnod,celvar,celeta,0.5/dy,0.5*pscale/dy,L);                        % P1
                        % Upper right cell
                        [L]=sigmayy_test(ky,ciB2,nody,celnod,celvar,celeta,0.5/dy,0.5*pscale/dy,L);                        % P1
                    end
                    % dSIGxy/dx
                    % x-distance between SIGxy-nodes
                    dx=abs(nodx(celnod(ci,4))-nodx(celnod(ci,2))); 
                    % Add SIGxy1
                    if(ciL2==0 && ciL1>0)
                        % vy between two subcells
                        %       |
                        %------ni1------
                        %       | ci
                        %  ciL1 +--vy3
                        %       | ciB1
                        %------ni2------
                        %       |
                        % Upper left node
                        [L,R]=sigmaxy_test(ky,celnod(ci,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,R,bupper,blower,bleft,bright,g);
                        % Lower left node
                        [L,R]=sigmaxy_test(ky,celnod(ciB1,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5/dx,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vy in standard position
                        %          |   ci 
                        %   ------ni1---vy3
                        %          |            
                        % Left node
                        [L,R]=sigmaxy_test(ky,celnod(ci,2),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-1/dx,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Add SIGxy2
                    if(ciR2==0 && ciR1>0)
                        % vy between two subcells
                        %         |
                        %  ------ni1------
                        %    ci   | 
                        %  --vy3--+ ciR1
                        %    ciB2 | 
                        %  ------ni2------
                        %         |
                        % Right top node
                        [L,R]=sigmaxy_test(ky,celnod(ci,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,R,bupper,blower,bleft,bright,g);
                        % Right bottom node
                        [L,R]=sigmaxy_test(ky,celnod(ciB2,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,0.5/dx,L,R,bupper,blower,bleft,bright,g);
                    else
                        % vy in standard position
                        %          |
                        %   vy3---ni1---
                        %     ci   |            
                        % Top node
                        [L,R]=sigmaxy_test(ky,celnod(ci,4),nodx,nody,nodcel,celnod,celvar,celres,nodeta,1/dx,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Right Part
%                     R(ky)=R(ky)-gy*(nodrho(celnod(ci,2))+nodrho(celnod(ci,4)))/2;
                    R(ky)=R(ky)-gy*celrho(ci,3);
                    % Mark processed equation
                    varyn(ky)=1;
                else
                    % dvy/dx stress gradient balance
                    dx=abs(nodx(celnod(ci,3))-nodx(celnod(ci,1)))/2; 
                    % Left subcel
                    if(ciB1>0)
%                         L(ky,ky)=-(nodeta(celnod(ciB1,1))+nodeta(celnod(ciB1,3)))/2/dx; % vy3
%                         L(ky,celvar(ciB1,2))=(nodeta(celnod(ciB1,1))+nodeta(celnod(ciB1,3)))/2/dx; % vy2
                        L(ky,ky)=-nodeta(celnod(ci,4))/dx; % vy3
                        L(ky,celvar(ciB1,2))=nodeta(celnod(ci,4))/dx; % vy2
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciB1,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciB1,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    % Right subcel
                    else
%                         L(ky,ky)=(nodeta(celnod(ciB2,1))+nodeta(celnod(ciB2,3)))/2/dx; % vy3
%                         L(ky,celvar(ciB2,2))=-(nodeta(celnod(ciB2,1))+nodeta(celnod(ciB2,3)))/2/dx; % vy2
                        L(ky,ky)=nodeta(celnod(ci,2))/dx; % vy3
                        L(ky,celvar(ciB2,2))=-nodeta(celnod(ci,2))/dx; % vy2
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciB2,1),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                        [L,R]=sigmaxy_vy_test(ky,celnod(ciB2,3),nodx,nody,nodcel,celnod,celvar,celres,nodeta,-0.5,L,R,bupper,blower,bleft,bright,g);
                    end
                    % Mark processed equation
                    varyn(ky)=1;
                end
            end
        end
        % End Equation for vy3------------------------
                

    end
end


% Solving the matrix
S=L\R;


% Computing stresses
exx=zeros(celnum,1);
eyy=zeros(celnum,1);
exy=zeros(nodnum,1);
sxx=zeros(celnum,1);
syy=zeros(celnum,1);
sxy=zeros(nodnum,1);
sxyyn=zeros(nodnum,1);
pr=zeros(nodnum,1);
% Going through all cells
for ci=1:1:celnum
    if(celdot(ci,1)==0)
        %
        %       |ciA1 |ciA2 | 
        %   ----+----vy2----+----
        %  ciL1 |  ci       | ciR1
        %   ---vx1   P5    vx4---
        %  ciL2 |           | ciR2
        %   ----+----vy3----+----
        %       |ciB1 |ciB2 |
        % Indexes of cells
        % Left
        ciL1=nodcel(celnod(ci,1),2); 
        ciL2=nodcel(celnod(ci,2),1);
        % Right
        ciR1=nodcel(celnod(ci,3),4);
        ciR2=nodcel(celnod(ci,4),3);
        % Above
        ciA1=nodcel(celnod(ci,1),3); 
        ciA2=nodcel(celnod(ci,3),1); 
        % Below
        ciB1=nodcel(celnod(ci,2),4);
        ciB2=nodcel(celnod(ci,4),2);
        
        % exx,eyy,sxx,syy,pr------------------------
        [exx1,sxx1,eyy1,syy1]=sigmaxxyy(ci,nodx,nody,celnod,celvar,celeta,S);
        exx(ci)=exx1;
        sxx(ci)=sxx1;
        eyy(ci)=eyy1;
        syy(ci)=syy1;
        pr(ci)=S(celvar(ci,5))*pscale;
        % End exx,eyy,sxx,syy,pr------------------------
        
        % Node 1 exy, sxy ------------------------
        ni=celnod(ci,1);
        if(sxyyn(ni)==0)
            if(nodx(ni)==g.xmin || nody(ni)==g.ymin || (ciL1>0 && ciA1>0))
                % Left/Upper boundary, not a middle node
                [exy1,sxy1]=sigmaxycalc(S,ni,nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=exy1;
                sxy(ni)=sxy1;
            else
                if(ciL1==0)
                % Left middle node
                % Upper node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciL2,3),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Lower node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciL2,4),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                else
                % Upper middle node
                % Left node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciA2,2),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Right node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciA2,4),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                end
            end
            % Mark current node
            sxyyn(ni)=1;
        end
       % End Node 1 exy, sxy ------------------------
        
        % Node 2 exy, sxy ------------------------
        ni=celnod(ci,2);
        if(sxyyn(ni)==0)
            if(nodx(ni)==g.xmin || nody(ni)==g.ymax || (ciL2>0 && ciB1>0))
                % Left/Upper boundary, not a middle node
                [exy1,sxy1]=sigmaxycalc(S,ni,nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=exy1;
                sxy(ni)=sxy1;
            else
                if(ciL2==0)
                % Left middle node
                % Upper node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciL1,3),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Lower node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciL1,4),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                else
                % Lower middle node
                % Left node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciB2,1),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Right node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciB2,3),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                end
            end
            % Mark current node
            sxyyn(ni)=1;
        end
       % End Node 2 exy, sxy ------------------------
        
        % Node 3 exy, sxy ------------------------
        ni=celnod(ci,3);
        if(sxyyn(ni)==0)
            if(nodx(ni)==g.xmax || nody(ni)==g.ymin || (ciR1>0 && ciA2>0))
                % Left/Upper boundary, not a middle node
                [exy1,sxy1]=sigmaxycalc(S,ni,nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=exy1;
                sxy(ni)=sxy1;
            else
                if(ciR1==0)
                % Right middle node
                % Upper node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciR2,1),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Lower node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciR2,2),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                else
                % Upper middle node
                % Left node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciA1,2),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Right node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciA1,4),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                end
            end
            % Mark current node
            sxyyn(ni)=1;
        end
       % End Node 3 exy, sxy ------------------------

        % Node 4 exy, sxy ------------------------
        ni=celnod(ci,4);
        if(sxyyn(ni)==0)
            if(nodx(ni)==g.xmax || nody(ni)==g.ymax || (ciR2>0 && ciB2>0))
                % Left/Upper boundary, not a middle node
                [exy1,sxy1]=sigmaxycalc(S,ni,nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=exy1;
                sxy(ni)=sxy1;
            else
                if(ciR1==0)
                % Right middle node
                % Upper node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciR1,1),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Lower node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciR1,2),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                else
                % Lower middle node
                % Left node
                [exy1,sxy1]=sigmaxycalc(S,celnod(ciB1,1),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                % Right node
                [exy2,sxy2]=sigmaxycalc(S,celnod(ciB1,1),nodx,nody,nodcel,nodeta,celnod,celvar,celres,bupper,blower,bleft,bright,g);
                exy(ni)=(exy1+exy2)/2;
                sxy(ni)=(sxy1+sxy2)/2;
                end
            end
            % Mark current node
            sxyyn(ni)=1;
        end
       % End Node 4 exy, sxy ------------------------
       
    end
end

% % Showing matrix structure
% figure(5)
% spy(L);
% 


end

