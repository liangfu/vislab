function sampling( dataFile, startNode, nPuppetsSamples )
% SAMPLING
% Visualize nPuppetsSamples samples from the DS model given a root (startNode) 
%
% Paper: "From Pictorial Structures to Deformable Structures", S. Zuffi, O.
% Freifeld, M. Black, CVPR 2012
% Authors: Silvia Zuffi, Department of Computer Science, Brown University
% Contact: zuffi@cs.brown.edu
% $Date: 2012 $
% $Revision: $
%
% Copyright 2011-2012, Brown University, Providence, RI. USA
%
%                          All Rights Reserved
%
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.
%
% For commercial uses contact the Technology Venture Office of Brown University
%
% The author and brown university disclaim all warranties with regard to
% this software, including all implied warranties of merchantability and
% fitness for any particular purpose.  In no event shall the author or
% brown university be liable for any special, indirect or consequential
% damages or any damages whatsoever resulting from loss of use, data or
% profits, whether in an action of contract, negligence or other tortious
% action, arising out of or in connection with the use or performance of
% this software.


SAVE = 0;

this = DS(dataFile);

fH = figure;

modelInd = this.modelForPartRef(startNode);

for i=1:nPuppetsSamples
  cla
  
  % Generate a random shape for the root, that we place in the middle of
  % the image with pi/2 angle
  mu_ab = mvnrnd(this.model{modelInd}.Mean, this.model{modelInd}.Covariance)';
  z_i = mu_ab(this.model{modelInd}.refCoeffInd);

  Bi = this.model{modelInd}.eigVectors{1};
  Sample = Bi * z_i + this.model{modelInd}.eigMean{1};
  X = Sample(1:2:end-2*this.nJointsPart(1));
  Y = Sample(2:2:end-2*this.nJointsPart(1));

  Oi0x = 0;
  Oi0y = 0;
  theta = pi/2;
  [ X Y ] = object_to_world(X, Y, theta, 0, 0);
  col = [ 0.0 0.7 0.0];
  plot(X, Y, 'Color', col, 'Linewidth', 2);
  hold on

  % Recursively generate children along the tree given the model
  get_and_plot_child(this, -1, startNode, z_i, theta, Oi0x, Oi0y, [], [], [], []);

  ylim([-0.5, 0.9]);
  axis equal
  axis ij

  drawnow
  if SAVE
      saveas(fH, sprintf('sample_%d', i), 'png');
  end
  pause(1)

end
end


function get_and_plot_child( this, parent, node, z_i, theta_i0, Oix_0, Oiy_0, z_i_s, theta_i0_s, Oix_0_s, Oiy_0_s)
  nSamples = 10;
  
  for a = 1:this.nParts
    if (this.jointConnectionMatrix(node, a) > 0)&&(a ~= parent)
      k = a;
      modelInd = find(ismember( this.modelName, this.modelIndex{node}{k} ) == 1);
      if any(modelInd)
          
        % Given the shape of the parent (z_i), compute the
        % conditioned model. The returned values are the mean of
        % the conditioned model. In the following, the Mean and
        % Covariance of the conditioned model are used to generate
        % samples
        [ z_j, sinTheta_ij, cosTheta_ij, Ojx, Ojy, Qjx, Qjy, Mean, Covariance ] = ...
            this.predict_part_conditioned( modelInd, z_i, [], [], [], 0 );

        Bj = this.model{modelInd}.eigVectors{2};
        M = this.model{modelInd}.eigMean{2};

        z_j_s = [];
        
        for s=1:nSamples
            try
                if any(z_i_s)
                    [ z_j_s(:,s), sinTheta_ij_s, cosTheta_ij_s, Ojx_s, Ojy_s ] =...
                        this.get_part_sample(this.model{modelInd}, Mean, Covariance, z_i_s(:,s));
                else
                    [ z_j_s(:,s), sinTheta_ij_s, cosTheta_ij_s, Ojx_s, Ojy_s ] =...
                        this.get_part_sample(this.model{modelInd}, Mean, Covariance, z_i);
                end
                Sample = Bj * z_j_s(:,s) + M;
                X = Sample(1:2:end-2*2);
                Y = Sample(2:2:end-2*2);
                
                theta_ji_s = atan2(sinTheta_ij_s, cosTheta_ij_s);
                if any(theta_i0_s)
                    theta_j0_s(s) = theta_ji_s + theta_i0_s(s);
                    [ Ojx_0_s(s) Ojy_0_s(s) ] = object_to_world(Ojx_s, Ojy_s, theta_i0_s(s), Oix_0_s(s), Oiy_0_s(s));
                else
                    theta_j0_s(s) = theta_ji_s + theta_i0;
                    [ Ojx_0_s(s) Ojy_0_s(s) ] = object_to_world(Ojx_s, Ojy_s, theta_i0, Oix_0, Oiy_0);
                end
                [ X Y ] = object_to_world(X, Y, theta_j0_s(s), Ojx_0_s(s), Ojy_0_s(s));

                col = [ 0.8 0.8 0.8];
                plot(X, Y, 'Color', col, 'Linewidth', 1);
                plot(X, Y, '--k', 'Linewidth', 1);
                hold on
            catch
                disp('UNABLE TO SAMPLE')
                Ojx_0_s = [];
                Ojy_0_s = [];
            end
        end
        
        % Plot the contour points
        Sample = Bj * z_j + M;
        X = Sample(1:2:end-2*2);
        Y = Sample(2:2:end-2*2);
        theta_ji = atan2(sinTheta_ij, cosTheta_ij);
        theta_j0 = theta_ji + theta_i0;
        [ Ojx_0 Ojy_0 ] = object_to_world(Ojx, Ojy, theta_i0, Oix_0, Oiy_0);
        [ X Y ] = object_to_world(X, Y, theta_j0, Ojx_0, Ojy_0);
        col = [ 0.0 0.0 0.0];
        plot(X(1:60), Y(1:60), 'Color', col, 'Linewidth', 2);
        hold on
        plot(X(61:end), Y(61:end), 'Color', col, 'Linewidth', 2);
        hold on

        [ Qjx_0 Qjy_0 ] = object_to_world(Qjx, Qjy, theta_i0, Oix_0, Oiy_0);
        X = Qjx_0;
        Y = Qjy_0;
        plot(X, Y, '*r');
        hold on
        
        get_and_plot_child( this, node, k, z_j, theta_j0, Ojx_0, Ojy_0, z_j_s, theta_j0_s, Ojx_0_s, Ojy_0_s );
      end
      
    end
  end
end


