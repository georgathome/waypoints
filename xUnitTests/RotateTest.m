classdef RotateTest < matlab.unittest.TestCase
    
    methods (Test)
        function testRotate(testCase)
            x = 1:10;
			y = zeros(size(x)) + 2;
			s = x - x(1);
			head = zeros(size(x));
			curv = zeros(size(x));
			
			obj = Waypoints.xy2Waypoints(x, y);
			
			dphi = pi/2;
			objR = rotate(obj, dphi);
			
			verifyEqual(testCase, objR.x, -y', 'AbsTol',1e-12);
			verifyEqual(testCase, objR.y, +x', 'AbsTol',1e-12);
			verifyEqual(testCase, objR.s, s');
			verifyEqual(testCase, objR.Head, head' + dphi);
			verifyEqual(testCase, objR.Curv, curv');
		end%fcn
		
    end
    
end%class