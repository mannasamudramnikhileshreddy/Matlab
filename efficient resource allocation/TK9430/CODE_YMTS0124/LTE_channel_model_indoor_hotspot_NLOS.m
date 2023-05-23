function Gain = LTE_channel_model_indoor_hotspot_NLOS( source_x,source_y, dest_x,dest_y)
     fc=2.6;
     sigma=8;
     distance_=sqrt((source_x-dest_x).^2+(source_y-dest_y).^2);
      rndnum=  sigma*randn;
      loss = 33.3*log10(distance_) + 11.5 + 40*log10(fc) +rndnum ;
      Gain = (1/(db2pow(loss)))*1e7;
end
