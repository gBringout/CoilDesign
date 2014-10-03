function [BX_a,BY_a,BZ_a] = Field2(BoucleCourant,I,X_Valeur,Y_Valeur,Z_Valeur)

% This file aim at calculating the B field accroding to biot savard equation 

%% Modification haute résolution
% % Z la position, en metre
% Z_Start =-0.15;
% Z_Stop = 0.15;
% Z_Step = 0.01;
% Z_NombreDeBoucle = (Z_Stop-Z_Start)/Z_Step + 1;
% 
% % X la position, en metre
% X_Start =-0.15;
% X_Stop = 0.15;
% X_Step = 0.01;
% X_NombreDeBoucle = (X_Stop-X_Start)/X_Step + 1;
% 
% % Y la position, en metre
% Y_Start = -0.15;
% Y_Stop = 0.15;
% Y_Step = 0.01;
% Y_NombreDeBoucle = (Y_Stop-Y_Start)/Y_Step + 1;
% 
% % Initialisation
% X_Valeur = zeros(1,X_NombreDeBoucle,'double');
% Y_Valeur = zeros(1,Y_NombreDeBoucle,'double');
% Z_Valeur = zeros(1,Z_NombreDeBoucle,'double');
% 
% for w=1:X_NombreDeBoucle
%     X_Valeur(1,w) = X_Start + (w-1)*X_Step;
% end
% for w=1:Y_NombreDeBoucle
%     Y_Valeur(1,w) = Y_Start + (w-1)*Y_Step;
% end
% for w=1:Z_NombreDeBoucle
%     Z_Valeur(1,w) = Z_Start + (w-1)*Z_Step;
% end

%% Calcule des constantes
X_NombreDeBoucle = size(X_Valeur,2);
Y_NombreDeBoucle = size(Y_Valeur,2);
Z_NombreDeBoucle = size(Z_Valeur,2);

NombreDeBoucle = size(BoucleCourant,2);
% DiametreDeCoil = a;
% I = 100;

%% Affichage des fils pour vérifier les données

% figure
% hold all;
% for w = 1:NombreDeBoucle
%     plot(BoucleCourant(w).Coord(1,:),BoucleCourant(w).Coord(2,:))
% end
% hold off;

%% Pré-allocation
BX_a = zeros(X_NombreDeBoucle,Y_NombreDeBoucle,Z_NombreDeBoucle);
BY_a = zeros(X_NombreDeBoucle,Y_NombreDeBoucle,Z_NombreDeBoucle);
BZ_a = zeros(X_NombreDeBoucle,Y_NombreDeBoucle,Z_NombreDeBoucle);

%% Constante
z=0;
totaltime = 0;
prev_lx = 0;
curr_lx = 0;
prev_ly = 0;
curr_ly = 0;
prev_lz = 0;
curr_lz = 0;
prev_coord_1u = 0;
coord_z = 0;
coord_theta = 0;
coord_1u = 0;
coord_2u = 0;
pi_180 = pi/180;
mu0 = 4*pi*10^-7;
mip4 = (mu0*I)/(4*pi);


    
%% Calcule du champs 3D
%han = waitbar(0,'Calcul en cours !');

disp('Calculation of the Field');
if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open
end
parfor x = 1:X_NombreDeBoucle;
% message = sprintf('%s%u%s%u%s%u%s%s%u%s %5.2f %s','Boucle X ',x,' sur ',X_NombreDeBoucle,'. Boucle Y (',Y_NombreDeBoucle,'). Boucle Z ','(',Z_NombreDeBoucle,'). Temps :',totaltime/60,'min');
% disp(message);
% tic;
% BoucleFaite = x;
% BoucleTotale = X_NombreDeBoucle;
% Avancement = BoucleFaite/BoucleTotale;
% waitbar(Avancement,han,sprintf('Avancement %2.1f %%. Temps restant : %5.2f min',Avancement*100, (totaltime/60) /Avancement - totaltime/60 ));
    tempX = zeros(Y_NombreDeBoucle,Z_NombreDeBoucle);
    tempY = zeros(Y_NombreDeBoucle,Z_NombreDeBoucle);
    tempZ = zeros(Y_NombreDeBoucle,Z_NombreDeBoucle);
    for y = 1:Y_NombreDeBoucle

        
        X_Val = X_Valeur(x);
        Y_Val = Y_Valeur(y);
        for z = 1:Z_NombreDeBoucle
            % Boucle pour passer dans tous le design de la coil
            Z_Val = Z_Valeur(z);
            somme_X = 0;
            somme_Y = 0;
            somme_Z = 0;
            
            for w=1:NombreDeBoucle
                sens = BoucleCourant(w).currentDirection;
                coord_z     = BoucleCourant(w).Coord(3,:);
                %coord_theta = BoucleCourant(w).Coord(2,:);
                coord_x = BoucleCourant(w).Coord(1,:);
                coord_y = BoucleCourant(w).Coord(2,:);
                
                prev_lz = coord_z(1);
                prev_lx = coord_x(1);
                prev_ly = coord_y(1);
%                 prev_lx = cos(coord_theta(1)*pi_180);
%                 prev_ly = sin(coord_theta(1)*pi_180);

                for u=2:size(BoucleCourant(w).Coord,2)

                    curr_lx = coord_x(u);
                    curr_ly = coord_y(u);
                    curr_lz = coord_z(u);
%                     curr_lx = cos(coord_theta(u)*pi_180);
%                     curr_ly = sin(coord_theta(u)*pi_180);
                    
                    lx = (curr_lx-prev_lx);
                    ly = (curr_ly-prev_ly);
                    lz = (curr_lz - prev_lz);
                    
                    rx = X_Val - ((curr_lx + prev_lx)/2);
                    ry = Y_Val - ((curr_ly + prev_ly)/2);
                    rz = Z_Val - ((curr_lz + prev_lz)/2);
                    
                    norm_vect = (sqrt(rx^2+ry^2+rz^2))^3;
                    coef_vect = sens*mip4 / norm_vect;
                    
                    somme_X = somme_X + coef_vect * (ly*rz-lz*ry);
                    somme_Y = somme_Y + coef_vect * (lz*rx-lx*rz);
                    somme_Z = somme_Z + coef_vect * (lx*ry-ly*rx);
                    
                    prev_lz      = curr_lz;
                    prev_lx      = curr_lx;
                    prev_ly      = curr_ly;
                end
            end
            tempX(y,z) = somme_X;
            tempY(y,z) = somme_Y;
            tempZ(y,z) = somme_Z;
        end
    end
    BX_a(x,:,:) = tempX;
    BY_a(x,:,:) = tempY;
    BZ_a(x,:,:) = tempZ;
%totaltime = totaltime + toc;
end

%% Sauvegarde des résultats
% 
% NomFichier = 'Test_DriveX_40x60elmts_Fields_.mat';
% disp(NomFichier);
% save(NomFichier,'X_Start','X_Stop','X_Step','X_Valeur','X_NombreDeBoucle',...
%                 'Y_Start','Y_Stop','Y_Step','Y_Valeur','Y_NombreDeBoucle',...
%                 'Z_Start','Z_Stop','Z_Step','Z_Valeur','Z_NombreDeBoucle',...
%                 'BoucleCourant','I','a',...
%                 'BX_a','BY_a','BZ_a');
