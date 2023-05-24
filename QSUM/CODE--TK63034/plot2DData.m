function plot2DData( fig, nSubplot,Base, data, centroid, minAxis, maxAxis, description )
    figure( fig );
    np = subplot( 1, 2, nSubplot );
    if centroid
        cla( np );
        hold on;
        plot(Base(:,:),'*b');
        plot( data( :, 1 ), data( :, 2 ), 'ro' );
        plot( centroid( :, 1 ), centroid( :, 2 ), 'kx', 'LineWidth', 3, 'MarkerSize', 15 );
        hold off;
    else
        hold on;
        plot(Base,'*b');
        plot( data( :, 1 ), data( :, 2 ), 'ro' );
        hold off;
    end
    if exist( 'description', 'var')
        title( description );
    end
    
    if exist( 'minAxis', 'var') && exist( 'maxAxis', 'var')
        width = 10;
        axis( [ minAxis( 1 ) + ( ( minAxis( 1 )*-1 )/abs( minAxis( 1 ) ) )*-width, ...
            maxAxis( 1 ) + ( ( maxAxis( 1 )*-1 )/abs( maxAxis( 1 ) ) )*-width, ...
            minAxis( 2 ) + ( ( minAxis( 2 )*-1 )/abs( minAxis( 2 ) ) )*-width, ...
            maxAxis( 2 ) + ( ( maxAxis( 2 )*-1 )/abs( maxAxis( 2 ) ) )*-width, ...
        ] )
    else
        axis square;
    end