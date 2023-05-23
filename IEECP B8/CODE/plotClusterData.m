function plotClusterData ( inputs, U, nC )
    colors = rand( nC, 3 );
    [ nInput, samples ] = size( inputs );
    figure( 'Name', 'Classified Data', 'NumberTitle', 'off' );
    for n = 1 : 1 : nInput
        [ ~, c ] = max( U( n, : ) );
        if samples == 3
            plot3( inputs( n, 1 ), inputs( n, 2 ), inputs( n, 3 ), 'Color', colors( c, : ), 'Marker', 'o' );
        else
            plot( inputs( n, 1 ), inputs( n, 2 ), 'Color', colors( c, : ), 'Marker', 'o' );
        end
        hold on;
    end
    if samples == 3
        grid on;
    end
    title( 'Classified Data' );