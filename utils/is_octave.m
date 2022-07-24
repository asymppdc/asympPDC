% function c = is_octave()
%
%  output: c = 0 if Matlab
%            > 0 if Octave environment

%subfunction that checks if we are in octave
function r = is_octave ()
  persistent x;
  if (isempty (x))
    x = exist ('OCTAVE_VERSION', 'builtin');
  end
  r = x;
end