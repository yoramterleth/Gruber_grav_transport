%% TIME_Headwall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% contribution of snow to glacier from surrounding terrain
%%% reads in climate variables, applies winstral, calculate gravitational 
%%% transport based on Gruber, 2007 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = GRUBER_headwall(H,C)

%% declare topography  

z = H.z ; 
%beta = H.beta ;
fnb = H.fnb ; 
Dmax = H.Dmax ; 

%% Load variables from previous timestep
Acc= zeros(size(H.Acc)) ; 
M = zeros(length(z(:,1)), length(z(1,:))); %H.M ;
Mstore = H.M ; 
Dstore = H.D ; 
D =  zeros(length(z(:,1)), length(z(1,:))); % H.D ;%
Fnb = zeros(length(z(:,1)), length(z(1,:)), 4) ;

%% load precipitation: 
I = H.P .* (H.T < C.rainsnowT-1);
I = I + H.P .* (C.rainsnowT-H.T+1)./2 .* ...
        (H.T < C.rainsnowT+1 & H.T > C.rainsnowT-1);

%% load the order of consideration 
ii = H.ii ;
jj = H.jj ; 

%% Is there new snow?
if nansum(I(:))> 0 %time.TCUR - min(IN.timelastsnow(:)) <= 1
    
    % print that there is grav. transport on this timestep
    disp('Snow Tranport...')

    %% gravitational transport     
    for i = 1:length(H.A(:,1))

       if ~isnan(H.A(i,3))       

        %% find location 
        ix = ii(i); % the x coord
        iy = jj(i); % the y coord

        %% calculate the inflow 
        M(iy,ix) = M(iy,ix) + I(iy,ix) ; % add the new precip to the mass 

        %% calculate outflow 

        % is snowmass larger than Dmax?
        if M(iy,ix) < Dmax(iy,ix)
           D(iy,ix) = M(iy,ix) ;                 
          else 
           D(iy,ix) = Dmax(iy,ix) ;
        end
        
        % register the relevant output: new deposits only
        % Acc(iy,ix) = D(iy,ix)-(Mstore(iy,ix)+I(iy,ix)); 
        
        for ik = 1:4

            if isnan(fnb(iy,ix,ik))
                Fnb(iy,ix,ik) = 0 ; 
            else 
            % calculate the ammount of mass to each surrounding cell
            Fnb(iy,ix,ik) = (M(iy,ix) - D(iy,ix)) * fnb(iy,ix,ik) ; 
            end 
        end
        
              
        if iy > min(H.jj(:)) && ix > min(H.ii(:))
        % adjust the mass of each surroudning cell
        M(iy-1,ix) = M(iy-1,ix) + Fnb(iy,ix,1) ; % inflow to cell above NB 1
        M(iy,ix-1) = M(iy,ix-1) + Fnb(iy,ix,2) ; % inflow to cell to the left NB 2
        M(iy,ix+1) = M(iy,ix+1) + Fnb(iy,ix,3) ; % inflow to cell to the right NB 3
        M(iy+1,ix) = M(iy+1,ix) + Fnb(iy,ix,4) ; % inflow to cell to the bottom NB 4

        % adjust the M follwoing outfow from current cell 
        M(iy,ix) = M(iy,ix) - sum(Fnb(iy,ix,:));

        % register the relevant output: new deposits only
        %Acc(iy,ix) = D(iy,ix)-(Mstore(iy,ix)+I(iy,ix)); 
        %Acc(iy,ix) = D(iy,ix)-I(iy,ix); 
        % register the relevant output: new deposits only
        Acc(iy,ix) = D(iy,ix)-(Mstore(iy,ix)+I(iy,ix));
        
        else 
            M(iy,ix) = 0 ;   % set M to zero at padded boundaries
            Acc(iy,ix) = 0 ; % set Acc to zero at padded boundaries
        end 
        
       
       end
    end


else 
    M = M + zeros(size(z));
    D = D + zeros(size(z));
    Acc = zeros(size(z)); 
    
end

%% Store output variables

H.D = D ;
H.M = M ; 
H.I = I ; 
H.Acc = Acc ; 
H.Acc(H.Acc<0)=0; % delete neg values produced when Dlim = 0 
H.Mout = M - (Mstore+I) ; % the change in M
H.Mout(H.Mout<0) = 0 ; 
H.Mstore = Mstore; 
H.Dstore = Dstore; 
end 
    


