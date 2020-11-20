function indFailed = test_scaleTangentToAxis(nbr)
	
	% define some test cases and their desired results
	tc = [...
		%solution	xLimits		yLimits		point		angle
		%
		% tangent point at lower left corner
		{[+10 0],	[0 10],		[0 10],		[0 0],		0};...
		{[+10 0],	[0 10],		[0 10],		[0 0],		pi/2};...
		{[0 -10],	[0 10],		[0 10],		[0 0],		pi};...
		{[0 -10],	[0 10],		[0 10],		[0 0],		3*pi/2};...
		%
		% tangent point at lower right corner
		{[0 -10],	[0 10],		[0 10],		[10 0],		0};...
		{[+10 0],	[0 10],		[0 10],		[10 0],		pi/2};...
		{[+10 0],	[0 10],		[0 10],		[10 0],		pi};...
		{[0 -10],	[0 10],		[0 10],		[10 0],		3*pi/2};...
		%
		% tangent point at upper right corner
		{[0 -10],	[0 10],		[0 10],		[10 10],	0};...
		{[0 -10],	[0 10],		[0 10],		[10 10],	pi/2};...
		{[+10 0],	[0 10],		[0 10],		[10 10],	pi};...
		{[+10 0],	[0 10],		[0 10],		[10 10],	3*pi/2};...
		%
		% tangent point at upper left corner
		{[+10 0],	[0 10],		[0 10],		[0 10],		0};...
		{[0 -10],	[0 10],		[0 10],		[0 10],		pi/2};...
		{[0 -10],	[0 10],		[0 10],		[0 10],		pi};...
		{[+10 0],	[0 10],		[0 10],		[0 10],		3*pi/2};...
		%
		% tangent point at upper left corner
		% k = tan(phi)
		% d = k*x - y(x) = 5 - 5*tan(phi)
		% r1 = sqrt(5^2 + (5-d)^2)
		% r2 = -r1
		{[sqrt(25+25*tan(pi/8)^2) -sqrt(25+25*tan(pi/8)^2)],...
					[0 10],		[0 10],		[5 5],		pi/8};
	];

% 			cell2struct(tc,{'sol','xLimits','yLimits','point','angle'},2);

	% handle input arguments
	if nargin < 1
		nbr = 1:size(tc,1);
	end%if

	% select test cases to perform
	tc = tc(nbr,:);
	nbrOfTestCases = size(tc,1);

	% pre-allocation
	isTestSuccesfull = false(nbrOfTestCases,1);

	% run test cases
	for i = 1:nbrOfTestCases

		[r1,r2] = Waypoints.scaleTangentToAxis(...
			tc{i,2},...
			tc{i,3},...
			tc{i,4},...
			tc{i,5});

		isTestSuccesfull(i) = all(tc{i,1} == [r1 r2]);

	end%for

	% evaluation and command line output
	nbrSuccesfull = sum(isTestSuccesfull);
	nbrFailed = sum(~isTestSuccesfull);
	fprintf('   Test results: \n')
	fprintf('     succesfull: %3d of %d.\n',nbrSuccesfull,nbrOfTestCases)
	fprintf('     failed:     %3d of %d.\n',nbrFailed,nbrOfTestCases)

	ind = 1:length(tc);
	indFailed = ind(~isTestSuccesfull);
	formatString = repmat('%d, ',1,nbrFailed);
	formatString = [formatString(1:end-2),'.'];
	if ~isempty(indFailed)
		fprintf(['     test cases failed: ',formatString,'\n'],...
			indFailed);
	end%if
	fprintf('\n');

end%fcn
