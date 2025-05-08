function shadeMPCregions(t, controllerMode,yl)
    % Adds shaded patches to indicate MPC-active regions on current plot
    wasMPC = false;

    for i = 1:length(controllerMode)
        if controllerMode(i) == 1 && ~wasMPC
            startIdx = i;
            wasMPC = true;
        elseif controllerMode(i) == 0 && wasMPC
            endIdx = i - 1;
            patch([t(startIdx), t(endIdx), t(endIdx), t(startIdx)], ...
                  [yl(1), yl(1), yl(2), yl(2)], ...
                  [0.8 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            wasMPC = false;
        end
    end

    % Handle case if still in MPC at the end
    if wasMPC
        endIdx = length(t);
        patch([t(startIdx), t(endIdx), t(endIdx), t(startIdx)], ...
              [yl(1), yl(1), yl(2), yl(2)], ...
              [0.8 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end
