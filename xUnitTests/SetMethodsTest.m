classdef SetMethodsTest < matlab.unittest.TestCase
    
    methods (Test)
        function testNonMonotonicS(testCase)
            x = 1:10;
			
			testCase.verifyError(@() Waypoints(x, x, fliplr(x), x, x), ...
				'WAYPOINTS:DecreasingS');
		end%fcn
		
    end
    
end%class