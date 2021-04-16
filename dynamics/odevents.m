function [value, isterminal, direction] = odevents(t, x)
% 
% Matt Werner (m.werner@vt.edu) - April 10, 2021
% 
% Mark special occasions throughout the flight profile as an event. This
% function is called by the Matlab ODE solver upon specifying the event
% function in the ODE settings.
% 
%    Inputs:
% 
%                 t - Integration time.
%                     Size: 1-by-1 (scalar)
%                     Units: s (seconds)
% 
%                 x - Integration state.
%                     Size: n-by-1 (vector)
%                     Units: ? (SI)
% 
%    Outputs:
% 
%             value - Determines whether an event is due according to its
%                     definition, which is (in general) dependent upon the
%                     earth, the rocket (its intended flight profile), the
%                     state, and time. 
%                     Note: An event occurs if any element crossing zero
%                           does so in the direction indicated by the
%                           corresponding element in 'direction'.
%                     Size: m-by-1 (vector)
%                     Units: ? (SI)
% 
%        isterminal - Decides whether or not the event should stop
%                     integration. Integration should be stopped at the end
%                     of an interval of ODE uniqueness. Such an end may
%                     correspond with any sudden discontinuity in the
%                     underlying dynamics determining the motion of the
%                     flight vehicle. For example, integration should be
%                     stopped upon staging the launch vehicle and at
%                     touchdown.
%                     Note: Stopping integration in the middle of an
%                           interval of uniqueness and continuing
%                           integration in the same interval does not
%                           affect uniqueness of that solution. Therefore,
%                           integration may be safely stopped and continued
%                           if uncertain that an action should halt
%                           integration.
%                     Size: m-by-1 (vector)
%                     Units: - (unitless)
% 
%         direction - Decides if a zero of 'value' is valid for marking an
%                     event. The zero's validity is determined by ensuring
%                     that it crosses zero in a particular specified
%                     direction. If the value crosses a zero but in the
%                     wrong direction, then there is no marking for an
%                     event. A direction of '0' marks events regardless of
%                     direction.
%                     Note: Scalars can cross 0 in two directions.
%                           1. (+1) - to + (upwards)
%                           2. (-1) + to - (downwards)
% 

% Mark the passage of every 8 minutes 40 seconds without stopping the
% integration
Tmark = 520;
value(1,1) = mod(t - Tmark, Tmark) - Tmark/2;
isterminal(1,1) = 0;
direction(1,1) = -1;