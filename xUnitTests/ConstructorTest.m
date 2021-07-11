classdef ConstructorTest < matlab.unittest.TestCase
    
    methods (Test)
        function testConstructorInputLength(testCase)
            x = 1:10;
			
			testCase.verifyError(@() Waypoints(x(1:end-1), x, x, x, x), ...
				'WAYPOINTS:Waypoints:numelInputs');
			
			testCase.verifyError(@() Waypoints(x, x(1:end-1), x, x, x), ...
				'WAYPOINTS:Waypoints:numelInputs');
			
			testCase.verifyError(@() Waypoints(x, x, x(1:end-1), x, x), ...
				'WAYPOINTS:Waypoints:numelInputs');
			
			testCase.verifyError(@() Waypoints(x, x, x, x(1:end-1), x), ...
				'WAYPOINTS:Waypoints:numelInputs');
			
			testCase.verifyError(@() Waypoints(x, x, x, x, x(1:end-1)), ...
				'WAYPOINTS:Waypoints:numelInputs');
		end%fcn
		
% 		function testConstructorInput
		
    end
    
end%class
