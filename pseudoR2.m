function pR2 = pseudoR2( realData, estimatedData, lambda )  
  EPS = 0.0000000001;
  nPoints1 = size( realData );
  nPoints2 = size( estimatedData );
  if ( nPoints1( 1 ) > nPoints1( 2 ) || nPoints2( 1 ) > nPoints2( 2 ) )
    error( 'Input data should be 1-dimensional vectors' );
  end
  if ( length( realData ) ~= length( estimatedData ) )
    error( 'Lengths of input real and estimated data are not the same' );
  else
    nPoints  = length( realData );
    if ( lambda == 0 )
      meanData = zeros( 1, nPoints ) + EPS;
    else
      meanData = zeros( 1, nPoints ) + lambda;
    end
  end
  estimatedData( estimatedData == 0 ) = EPS;
  division1 = realData./estimatedData;
  division2 = realData./meanData;
  division1( division1 == 0 ) = EPS;
  division2( division2 == 0 ) = EPS;
  sum1 = sum( realData.*log( division1 ) - ( realData - estimatedData ) );
  sum2 = sum( realData.*log( division2 ) - ( realData - meanData ) );
  if ( sum2 == 0 )
    sum2 = EPS;
  end
  pR2 = 1 - sum1/sum2;