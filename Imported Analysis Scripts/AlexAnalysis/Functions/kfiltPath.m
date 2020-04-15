function kfiltPath(paths)
    for p = paths'
        s = load(p{1});
        pos = s.pos.uninterp;
        kpos = nan(size(pos));
        
        param = getDefaultParameters();
%         param.motionModel = 'ConstantVelocity';
%         param.initialEstimateError = param.initialEstimateError(1:2);
%         param.motionNoise          = param.motionNoise(1:2);
        
        kalmanFilter = configureKalmanFilter(param.motionModel, ...
          pos(:,1)', param.initialEstimateError, ...
          param.motionNoise, param.measurementNoise);
        for i = 1:length(pos(1,:))
            if ~isnan(pos(1,i))
                if i~=1
                    predict(kalmanFilter);
                end
                kpos(:,i) = correct(kalmanFilter, pos(:,i)');
            else
                kpos(:,i) = predict(kalmanFilter);
            end
        end
    end
end

function param = getDefaultParameters
  param.motionModel           = 'ConstantAcceleration';
  param.initialLocation       = 'Same as first detection';
  param.initialEstimateError  = 1E5 * ones(1, 3);
  param.motionNoise           = [25, 10, 1];
  param.measurementNoise      = 25;
  param.segmentationThreshold = 0.05;
end
