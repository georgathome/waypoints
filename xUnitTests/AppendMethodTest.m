classdef AppendMethodTest < matlab.unittest.TestCase
    
    methods (Test)
        function testNonMonotonicS(testCase)
            x = 1:10;
			obj = Waypoints.xy2Waypoints(x, x);
			objA = obj.append(obj);
			
			testCase.verifyTrue(~any(diff(objA.s < 0)))
		end%fcn
		
    end
    
end%class