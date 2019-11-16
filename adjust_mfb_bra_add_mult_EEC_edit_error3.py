import matplotlib as mpl
from datetime import datetime,timedelta
from mpl_toolkits.basemap import Basemap 
from statistics import mean
import matplotlib.pyplot as plt
import wradlib as wrl
import numpy as np
import pandas as pd
import warnings,os,sys,math,openpyxl,csv
from collections import Counter
from scipy.interpolate import griddata,Rbf
from datetime import datetime,timedelta
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=RuntimeWarning)

#   melihat errornya di line brapa
def errlog():
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno, exc_obj)
    pass

#   
def radar_time(f):
    try:time=datetime(int(f[3:7]),int(f[7:9]),int(f[9:11]),int(f[11:13]),int(f[13:15]))
    except:
        try:time=datetime(int(f[8:12]),int(f[12:14]),int(f[14:16]),int(f[17:19]),int(f[19:21]))
        except:
            try:time=datetime(int(f[4:8]),int(f[8:10]),int(f[10:12]),int(f[13:15]),int(f[15:17]))
            except:
                time=datetime(int(f[0:4]),int(f[4:6]),int(f[6:8]),int(f[8:10]),int(f[10:12])) 
    return time


def collectdata(datapath,timeacc,vendor):
    '''Collect data in the time window accumulation'''
    print ('\n----------------------------------------------------------------------')
    if vendor=='GEMA':strtimeacc=timeacc.strftime("%Y%m%d%H")
    elif vendor=='EEC':strtimeacc=timeacc.strftime("%Y%m%d-%H")
    print ('Collect base time : {}*'.format(strtimeacc))
    listpcdfpath,times,sfpath,stime=list(),list(),list(),list()
    for filelist in os.listdir(datapath):
        if strtimeacc in filelist:
            fpath=datapath+'/'+filelist
            listpcdfpath.append(fpath)
    '''Defining time in each file for calculating qpe'''
    for i in listpcdfpath:
        ff=i[(len(datapath)+1):]
        t=radar_time(ff)
        times.append(t)
    sfpath=sorted(listpcdfpath)
    stime=sorted(times)
    return sfpath,stime


#   
def sizegrid(fpath,res,rmax):
    f=wrl.util.get_wradlib_data_file(fpath)
    raw=wrl.io.read_generic_netcdf(f)

    llon = float(raw['variables']['longitude']['data'])
    llat = float(raw['variables']['latitude']['data'])
    lalt = float(raw['variables']['altitude']['data'])
        
    sitecoords=(llon,llat,lalt)
    res_coords=res/111229.
    xmax,xmin=llon+(rmax/111229),llon-(rmax/111229)
    ymax,ymin=llat+(rmax/111229),llat-(rmax/111229)
    n_grid=np.floor(((xmax-xmin)/res_coords)+1)
    x_grid=np.linspace(xmax,xmin,n_grid)
    y_grid=np.linspace(ymax,ymin,n_grid) 
    grid_xy=np.meshgrid(x_grid,y_grid)
    size_grid=grid_xy[0]
    return sitecoords,size_grid,x_grid,y_grid

#   
def EEC(fpath,res,nelevation,rmax):
    f=wrl.util.get_wradlib_data_file(fpath)
    raw=wrl.io.read_generic_netcdf(f)

    llon = float(raw['variables']['longitude']['data'])
    llat = float(raw['variables']['latitude']['data'])
    lalt = float(raw['variables']['altitude']['data'])

    sitecoords = (llon, llat,lalt)
    res_coords=res/111229.
    
    xmax,xmin=llon+(rmax/111229),llon-(rmax/111229)
    ymax,ymin=llat+(rmax/111229),llat-(rmax/111229)
    n_grid=np.floor(((xmax-xmin)/res_coords)+1)
    x_grid=np.linspace(xmax,xmin,n_grid)
    y_grid=np.linspace(ymax,ymin,n_grid)

    all_data = np.zeros((len(x_grid),len(y_grid)))

    #load range data
    r_all = raw['variables']['range']['data'] 
    r=r_all

    #load elevation
    nelevangle = np.size(raw['variables']['fixed_angle']['data'])
    sweep_start_idx = raw['variables']['sweep_start_ray_index']['data']
    sweep_end_idx = raw['variables']['sweep_end_ray_index']['data']
    n_azi = np.size(raw['variables']['azimuth']['data'])/nelevangle
    
    
    #flag option for loading EEC radar data
    sweep_start_idx = raw['variables']['sweep_start_ray_index']['data']
    sweep_end_idx = raw['variables']['sweep_end_ray_index']['data']
    try:
        if raw['gates_vary']=='true':
            ray_n_gates=raw['variables']['ray_n_gates']['data']
            ray_start_index=raw['variables']['ray_start_index']['data']
            flag='true'
        elif raw['gates_vary']=='false':
            flag='false'
    except :
        if raw['n_gates_vary']=='true':
            ray_n_gates=raw['variables']['ray_n_gates']['data']
            ray_start_index=raw['variables']['ray_start_index']['data']
            flag='true'
        elif raw['n_gates_vary']=='false':
            flag='false'
            
    if nelevation>nelevangle:
        nelevation=nelevangle

##    # Ekstrak v data untuk interferensi   
##    if flag == 'false':
##        data_v = raw['variables']['VELH']['data'][sweep_start_idx[0]:sweep_end_idx[0], :]    
##    else: #flag = true
##        data_v = np.array([])
##        n_azi = sweep_end_idx[0]-sweep_start_idx[0]        
##        for ll in range(sweep_start_idx[0],sweep_end_idx[0]):
##            data_v = np.append(data_v,raw['variables']['VELH']['data'][ray_start_index[ll]:ray_start_index[ll+1]])
##        data_v = data_v.reshape((n_azi,ray_n_gates[sweep_start_idx[0]]))
##    datav=data_v
##    datav[(datav==-32768.)]=-999.0


    for i in range(nelevation):
        elevation=float('{0:.1f}'.format(raw['variables']['fixed_angle']['data'][i]))        
        strsweep=str(i+1) 
        print ('Extracting radar data : SWEEP-'+strsweep +' at Elevation Angle '+str(elevation)+'  deg ...')

        # Load azimuth data  
        azi = raw['variables']['azimuth']['data'][sweep_start_idx[i]:sweep_end_idx[i]]   
            
        # Load radar data     
        if flag == 'false':
            data_ = raw['variables']['DBZH']['data'][sweep_start_idx[i]:sweep_end_idx[i], :]
            # create range array
            r = r_all    
        else: #flag = true                
            data_ = np.array([])
            n_azi = sweep_end_idx[i]-sweep_start_idx[i]        
            try:
                for ll in range(sweep_start_idx[i],sweep_end_idx[i]):
                    data_ = np.append(data_,raw['variables']['DBZH']['data'][ray_start_index[ll]:ray_start_index[ll+1]])
                data_ = data_.reshape((n_azi,ray_n_gates[sweep_start_idx[i]]))
            except:
                for ll in range(sweep_start_idx[i],sweep_end_idx[i]):
                    data_ = np.append(data_,raw['variables']['UH']['data'][ray_start_index[ll]:ray_start_index[ll+1]])
                data_ = data_.reshape((n_azi,ray_n_gates[sweep_start_idx[i]]))
            # create range array
            r = r_all[0:ray_n_gates[sweep_start_idx[i]]]

        data=data_
        clutter=wrl.clutter.filter_gabella(data, tr1=6, n_p=6, tr2=1.3, rm_nans=False)
        data_noclutter=wrl.ipol.interpolate_polar(data, clutter, Interpolator = wrl.ipol.Linear)
        data=data_noclutter
        

        # Ekstrak v data untuk interferensi
        filter=[0,1]
##        filter=[5]
        if i in filter:
            print ('Activating interference filter in elevation {}'.format(elevation))
            if flag == 'false':
                data_v = raw['variables']['VELH']['data'][sweep_start_idx[i]:sweep_end_idx[i], :]    
            else: #flag = true
                data_v = np.array([])
                n_azi = sweep_end_idx[i]-sweep_start_idx[i]        
                for ll in range(sweep_start_idx[i],sweep_end_idx[i]):
                    data_v = np.append(data_v,raw['variables']['VELH']['data'][ray_start_index[ll]:ray_start_index[ll+1]])
                data_v = data_v.reshape((n_azi,ray_n_gates[sweep_start_idx[i]]))
            datav=data_v
            datav[(datav==-32768.)]=-999.0
            data[datav==-999.0]=np.nan


        polargrid=np.meshgrid(r,azi)
        x,y,z=wrl.georef.polar2lonlatalt_n(polargrid[0],polargrid[1],elevation,sitecoords)

        grid_xy=np.meshgrid(x_grid,y_grid)
        xgrid=grid_xy[0]
        ygrid=grid_xy[1]
        grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()
        xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
        radius=r[np.size(r)-1]
        center=[x.mean(),y.mean()]
        # Interpolate polar coordinate data to cartesian cordinate
        # Option : Idw, Linear, Nearest
        gridded = wrl.comp.togrid(xy, grid_xy, radius, center, data.ravel(), wrl.ipol.Linear)
        gridded_data = np.ma.masked_invalid(gridded).reshape((len(x_grid), len(y_grid)))
        all_data=np.dstack((all_data,gridded_data))

    data=np.nanmax(all_data[:,:,:],axis=2)
    data[data<0]=np.nan;data[data>100]=np.nan
    xymaxmin=(xmax,xmin,ymax,ymin)    
    return data,sitecoords,xgrid,ygrid,xymaxmin
    
def calculate_qpe(datapath,sfpath,stime,res,nelevation,size_grid,rmax):
    qpe_data=np.zeros((np.shape(size_grid)[0],np.shape(size_grid)[1]))
    for file in sfpath:
        print ('Extract data : {}'.format(file[len(datapath)+1:]))
        dbz,sitecoords,xgrid,ygrid,xymaxmin=EEC(file,res,nelevation,rmax)
        z=wrl.trafo.idecibel(dbz)
        r=wrl.zr.z2r(z,a=250,b=1.2)
        qpe_data=np.dstack((qpe_data,r))
    qpe_data[qpe_data==0.036]=np.nan
    qpe_data=np.delete(qpe_data,0,axis=2)
    print ('\nCalculating qpe...')
    x,y,z=np.shape(qpe_data)
    qpe,A=np.zeros((x,y)),np.zeros((x,y))
    sec_in_hours=timedelta(hours=1).total_seconds()
    delta=list()
    for i in range(len(stime)):
        if i != 0:
            d=stime[i]-stime[i-1]
            ts=d.total_seconds()
            dd=ts/sec_in_hours
            delta.append(dd)
    for k in range(z):
        if k != 0:
            for i in range(x):
                for j in range(y):               
                    A[i,j]=delta[k-1]*((qpe_data[i,j,k]+qpe_data[i,j,(k-1)])/2)
            qpe=np.dstack((qpe,A))

    
    qpe=np.delete(qpe,0,axis=2)
    data_qpe=np.sum(qpe,axis=2)
    
##    mask=Counter(data_qpe.flatten()).most_common(1);value,count=mask[0][0],mask[0][1]
##    data_qpe[data_qpe<0.1]=np.nan
##    if (int(count)>100):data_qpe[data_qpe==value]=np.nan
    data_qpe[data_qpe==data_qpe[0,0]]=np.nan
    print ('Finished calculating qpe1h.\n')
    return data_qpe


def extract_point(lat_arg,lon_arg,lat_radar,lon_radar,qpe):
    array_lat=np.asarray(lat_radar)
    array_lon=np.asarray(lon_radar)
    idxlat=(np.abs(array_lat-lat_arg)).argmin()
    idxlon=(np.abs(array_lon-lon_arg)).argmin()
    radar_point=qpe[idxlon,idxlon]
    point_lat=lat_radar[idxlat]
    point_lon=lon_radar[idxlon]

    return point_lat,point_lon,radar_point

##  gauge yang kurang dari 80%
def read_gauge_not_use(filegaugenotuse):
    gauge=csv.reader(open(filegaugenotuse),delimiter='\n')
    gauge_not_use_opr_blck=[]
    for row in gauge:
        gauge_not_use_opr_blck.append(row[0])
    return gauge_not_use_opr_blck

##  gauge yang terpakai dalam range radar
def read_gauge_use(filemetadata,provinsi,gauge_not_use,sitecoords,rmax_use):
    llon,llat,lalt=sitecoords
    xmax,xmin=llon+(rmax_use/111229),llon-(rmax_use/111229)
    ymax,ymin=llat+(rmax_use/111229),llat-(rmax_use/111229)
    gauge_use=[]
    metadata=csv.reader(open(filemetadata),delimiter=',')
    for row in metadata:
        try:
            lon_arg=float(row[4]);lat_arg=float(row[3])
            if row[2] in provinsi:
                if row[0] in gauge_not_use:pass
                elif ymin < lat_arg < ymax and \
                    xmin < lon_arg < xmax :
                    gauge_use.append(row[0])
        except:pass
    gauge_use = [x for x in gauge_use if x != '']
    return gauge_use

    print()
    
def datapoint(filemetadata,fileaws,provinsi,gauge_use,lat_radar,lon_radar,qpe,timeacc,data_unrealible):
    wb = openpyxl.load_workbook(fileaws)    
    gauge_acc,radar_acc=np.array([]),np.array([])
    gauge_lon,gauge_lat=np.array([]),np.array([])
    metadata=csv.reader(open(filemetadata),delimiter=',')
    for i in range(len(provinsi)):
        for j in range(len(gauge_use)):
            ws=wb[provinsi[i]]
            rows=ws.max_row;columns=ws.max_column
            for k in range(columns):
                if gauge_use[j]==str(ws.cell(row=1,column=k+1).value):
                    index_column=k+1
                    init=0.
                    for l in range(rows):
                        try:
                            if timeacc.strftime("%Y-%m-%d %H") in str(ws.cell(row=l+1,column=1).value):
                                value=float(str(ws.cell(row=l+1,column=index_column).value))
                                if value>data_unrealible : value=0.
                                init=init+value
                        except:pass
                    gauge_acc=np.append(gauge_acc,init)    
    for row in metadata:
        for j in range(len(gauge_use)):
            if gauge_use[j]==row[0]:
                lon_arg=float(row[4]);lat_arg=float(row[3])
                gauge_lon=np.append(gauge_lon,lon_arg)
                gauge_lat=np.append(gauge_lat,lat_arg)
                point_lat,point_lon,radar_point=extract_point(lat_arg,lon_arg,lat_radar,lon_radar,qpe)
                radar_acc=np.append(radar_acc,radar_point)
    return gauge_acc,radar_acc,gauge_lon,gauge_lat
#    print (gauge_acc,radar_acc,gauge_lon,gauge_lat)
##    df = pd.DataFrame(radar_acc)
##    df.to_csv('radar_acc_2019032613.csv', sep=',', index=False)
    
def calculate_distance(pointA,pointB):
    # approximate radius of earth in km
    R = 6373.0
    lat1 = math.radians(pointA[1])
    lon1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[1])
    lon2 = math.radians(pointB[0])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    return distance

##    print 'sampai disini'    
def adjust(gauge_acc,radar_acc,gauge_lon,gauge_lat,sitecoords,rmax_use,x_grid,y_grid,interpmethod,maxthreshold):
    '''Mean Field Bias'''
    print ('Calculating Mean Field Bias correction factor')
    x_target,y_target=np.meshgrid(x_grid,y_grid)
##    print 'sampai disini'
    correction_mfb=np.array([])
#    print 'sampai disini'
    for i in range(len(gauge_acc)):
        if str(gauge_acc[i])=='nan' or str(gauge_acc[i])=='inf':gauge_acc[i]=0.
        if str(radar_acc[i])=='nan' or str(radar_acc[i])=='inf':radar_acc[i]=0.
        c=gauge_acc[i]/radar_acc[i]
        if str(c)=='nan' or str(c)=='inf':c=1
##        print ('sampai disini')
        correction_mfb=np.append(correction_mfb,c)
##        print 'sampai disini'
    correction_mfb[correction_mfb<0.2]=1.
##    print 'sampai disini'
    correction_mfb[correction_mfb>maxthreshold]=maxthreshold
##    print 'sampai disini'
    
    '''Brandes Spatial Adjustment'''
    print ('Calculating Brandes Spatial Adjustment correction factor')
    llon,llat,lalt=sitecoords
    number_of_gauges=float(len(gauge_acc))
    miu=number_of_gauges/(np.pi*pow(rmax_use/1000.,2))
    k=pow((2*miu),-1.)
    weight=np.array([]);correction_bra=np.array([])
    for i in range(len(gauge_acc)):
        pointA=(llon,llat);pointB=(gauge_lon[i],gauge_lat[i])
        d=calculate_distance(pointA,pointB)
        wi=math.exp(-(pow(d,2)/k))
        weight=np.append(weight,wi)
    sum_weight=np.nansum(weight)
    for i in range(len(gauge_acc)):
        c_bra=(weight[i]*(gauge_acc[i]/radar_acc[i]))/sum_weight
        if str(c_bra)=='nan' or str(c_bra)=='inf':c_bra=1
        correction_bra=np.append(correction_bra,c_bra)
    correction_bra[correction_bra<0.2]=1.
#    correction_bra[correction_bra>maxthreshold]=maxthreshold

    '''Additive Correction'''
    correction_add=np.array([])
    for i in range(len(gauge_acc)):
        c_add=gauge_acc[i]-(radar_acc[i]*correction_mfb[i])
        correction_add=np.append(correction_add,c_add)
        
    '''Interpolate point correction to grid'''
    if interpmethod=='griddata':
        print ('\nInterpolate correction factor into grid (scipy.interp.griddata)')
        grid_mfb=griddata((gauge_lon,gauge_lat),correction_mfb,(x_target,y_target),method='cubic')
        grid_bra=griddata((gauge_lon,gauge_lat),correction_bra,(x_target,y_target),method='cubic')
        grid_add=griddata((gauge_lon,gauge_lat),correction_add,(x_target,y_target),method='cubic')
        grid_gauge=griddata((gauge_lon,gauge_lat),gauge_acc,(x_target,y_target),method='cubic')
    elif interpmethod=='Rbf':
        print ('\nInterpolate correction factor into grid (scipy.interp.Rbf)')
        rbf_gauge = Rbf(gauge_lon, gauge_lat, gauge_acc, function='linear')
        rbf_mfb = Rbf(gauge_lon, gauge_lat, correction_mfb, function='linear')
        rbf_bra = Rbf(gauge_lon, gauge_lat, correction_bra, function='linear')
        rbf_add = Rbf(gauge_lon, gauge_lat, correction_add, function='linear')
        grid_gauge = rbf_gauge(x_target, y_target)
        grid_mfb = rbf_mfb(x_target, y_target)
        grid_bra = rbf_bra(x_target, y_target)
        grid_add = rbf_add(x_target, y_target)
        
    grid_bra=np.nan_to_num(grid_bra);grid_mfb=np.nan_to_num(grid_mfb)
    grid_bra[grid_bra==0.]=1;grid_mfb[grid_mfb==0.]=1.
    grid_bra[grid_bra<0.2]=1;grid_mfb[grid_mfb<0.2]=1.
    return correction_mfb,correction_bra,correction_add,grid_mfb,grid_bra,grid_add,grid_gauge

def verification(gauge_acc,radar_acc,corr,adjtype):
    init_mae=0;init_me=0;init_rmse=0.
    init_mae_adj=0;init_me_adj=0;init_rmse_adj=0.
    bias=np.array([])
    for i in range(len(gauge_acc)):
        if str(gauge_acc[i])=='nan' or str(gauge_acc[i])=='inf':gauge_acc[i]=0.
        if str(radar_acc[i])=='nan' or str(radar_acc[i])=='nan':radar_acc[i]=0.
        
        err=radar_acc[i]-gauge_acc[i]
        abs_err=abs(radar_acc[i]-gauge_acc[i])
        bias=np.append(bias,abs_err)
        err_pow2=(radar_acc[i]-gauge_acc[i])**2

        if adjtype=='multiplicative':
            abs_err_adj=abs((radar_acc[i]*corr[i])-gauge_acc[i])
            err_adj=(radar_acc[i]*corr[i])-gauge_acc[i]
            err_pow2_adj=((radar_acc[i]*corr[i])-gauge_acc[i])**2
        elif adjtype=='additive':
            abs_err_adj=abs((radar_acc[i]+corr[i])-gauge_acc[i])
            err_adj=(radar_acc[i]+corr[i])-gauge_acc[i]
            err_pow2_adj=((radar_acc[i]+corr[i])-gauge_acc[i])**2
            
        init_mae=init_mae+abs_err
        init_mae_adj=init_mae_adj+abs_err_adj
        init_me=init_me+err
        init_me_adj=init_me_adj+err_adj
        init_rmse=init_rmse+err_pow2
        init_rmse_adj=init_rmse_adj+err_pow2_adj
    mae=init_mae/float(len(gauge_acc))
    mae_adj=init_mae_adj/float(len(gauge_acc))
    me=init_me/float(len(gauge_acc))
    me_adj=init_me_adj/float(len(gauge_acc))
    rmse=math.sqrt((init_rmse/float(len(gauge_acc))))
    rmse_adj=math.sqrt((init_rmse_adj/float(len(gauge_acc))))
    print ('\nMAE : {0:.4f}'.format(mae))
    print ('ME : {0:.4f}'.format(me))
    print ('RMSE : {0:.4f}'.format(rmse))
    print ('MAE adj : {0:.4f}'.format(mae_adj))
    print ('ME adj  : {0:.4f}'.format(me_adj))
    print ('RMSEadj : {0:.4f}'.format(rmse_adj))
    return mae,me,rmse,mae_adj,me_adj,rmse_adj,bias

def MAKE_CMAP(colors, position=None, bit=False):   
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
        
    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap
##    print 'sampai disini'
    
def plotonmap(lat,lon,shppath,data,title,fileimage,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type): 
    '''Masking'''
    mask_data_=data
    mask_ind=np.where(mask_data_ <= np.nanmin(mask_data_))
    mask_data_[mask_ind]=np.nan
    data=np.ma.array(mask_data_,mask=np.isnan(mask_data_))
    
    center_lon,center_lat=lon.mean(),lat.mean()
    left_lon,right_lon=lon.min(),lon.max()
    top_lat,down_lat=lat.max(),lat.min()
    m=Basemap(llcrnrlat=down_lat,urcrnrlat=top_lat,\
                llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='i')

    x,y=np.meshgrid(lon,lat)
    x1,y1=m(x,y)
    x0,y0=m(center_lon,center_lat) ##
    
    if type=='rain':
        Rainbow_Rainrate=[(0,0,199),(0,121,255),(50,200,255),(120,235,255),\
                        (255,255,255),(255,247,192),(255,229,0),(255,115,0),\
                        (255,63,0),(200,0,0),(150,0,0),(110,0,0)] 
        clevs_R = [0.1,1.0,2.0,5.0,7.0,9.0,10,12,15,20,30,50]
##        clevs_R = [0.1,1.0,5.0,9.0,15,20,25,30,40,50,60,100]
        color_rain=Rainbow_Rainrate
        mycmap_rain=MAKE_CMAP(color_rain,bit=True)   
        plt.figure(figsize=(10,8))
        for i in range(len(gauge_acc)):
            xpt,ypt=m(gauge_lon[i],gauge_lat[i])
            if title[:4]=='Gaug':
                m.plot(xpt,ypt,'bo',markersize=1)
                plt.text(xpt,ypt,gauge_acc[i],fontsize=5,color='black')      
        m.contourf(x1, y1, data,clevs_R,colors=Rainbow_Rainrate,extend='max')
        m.colorbar(ticks=clevs_R,label='Rain Accumulation (mm)')
    elif type=='correction':
        plt.figure(figsize=(10,8))
        for i in range(len(gauge_acc)):
            xpt,ypt=m(gauge_lon[i],gauge_lat[i])
            m.plot(xpt,ypt,'bo',markersize=1)
            if title[:4]=='Mean':
                plt.text(xpt,ypt,correction_mfb[i],fontsize=5,color='black')
            elif title[:4]=='Bran':
                plt.text(xpt,ypt,correction_bra[i],fontsize=5,color='black')
            
        m.contourf(x1, y1, data,cmap='rainbow')
        m.colorbar(label='Correction')
    m.plot(x0,y0,'ko',markersize=3) ##
    m.drawparallels(np.arange(math.floor(y.min()),math.ceil(y.max()),1.5),labels=[1,0,0,0],linewidth=0.2,fontsize=5,color='black')
    m.drawmeridians(np.arange(math.floor(x.min()),math.ceil(x.max()),1.5),labels=[0,0,0,1],linewidth=0.2,fontsize=5,color='black')  
    #m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2500, verbose= True)
    m.drawcoastlines()
##    m.arcgisimage(service='World_Shaded_Relief', xpixels = 1500, verbose= True)
    m.readshapefile(shppath,'Kabupaten',linewidth=0.4,color='gray')
    plt.title(title,weight='bold',fontsize=15)
    plt.savefig(fileimage,bbox_inches='tight',dpi=200,pad_inches=0.5)
    plt.close()
    #plt.show()
    #m.bluemarble()

def write_csv(gauge_use,gauge_acc,radar_acc,timeacc,csvpath,correction_mfb,correction_bra,correction_add):
    csv_file='{}/{}'.format(csvpath,timeacc.strftime("%Y%m%d%H00.csv"))
    with open (csv_file,'w') as f:
        for i in range(len(gauge_use)):
            f.write('{0},{1:.4f},{2:.4f},{3:.4f},{4:.4f},{5:.4f}\n'.format(gauge_use[i],gauge_acc[i],radar_acc[i],\
                                                                           np.nan_to_num(radar_acc[i])*correction_mfb[i],\
                                                                           np.nan_to_num(radar_acc[i])*correction_bra[i],\
                                                                           np.nan_to_num(radar_acc[i])+correction_add[i]))
        f.close()

'''Main Program'''
##os.environ['WRADLIB_DATA']='D:/wradlib-data-master'
os.environ['WRADLIB_DATA']='D:/6.RADAR/handsontraining/RADAR/wradlib-data-master'
time_initial=datetime(2018,12,12,9)
time_window=1
for i in range(time_window):
    timeacc=time_initial+timedelta(hours=i)
##    csvpath='D:/BISMILLAH/proses/csv'
    csvpath='D:/BISMILLAH/proses/adjustradar/data/JAK/20181212-09'
##    radarpath='D:/BISMILLAH/proses/adjustradar/data/JAK/maret'
    radarpath='D:/BISMILLAH/proses/adjustradar/data/JAK/20181212-09'
##    imagepath='D:/BISMILLAH/proses/adjustradar/output'
##    imagepath='D:/adjustradar/output/20190326_rbf_linear'
    imagepath='D:/BISMILLAH/proses/adjustradar/data/JAK/20181212-09'
    shppath='D:/BISMILLAH/proses/adjustradar/shp_lite/Kabupaten'
    gaugepath='D:/BISMILLAH/proses/adjustradar/excel/coba_error/edit'
##    print 'sampai disini'
    filemetadata='D:/BISMILLAH/proses/adjustradar/metadata.txt'
##    filegaugenotuse='G:/adjustradar/gauge_not_use_JAK_20190326.txt'
    fileaws='{}/rr_indonesia_{}.xlsx'.format(gaugepath,timeacc.strftime("%Y%m%d"))
    fileimage_cont='{}/contingency{}_qpe1h.png'.format(gaugepath,timeacc.strftime("%Y%m%d%H00"))
    fileimage_scatter='{}/scatter-radar-gauge{}_qpe1h.png'.format(imagepath,timeacc.strftime("%Y%m%d%H00"))
    provinsi=['DKI Jakarta','Jawa Barat','Banten'];vendor='EEC';site='JAK';nelevation=4;res=1000.;maxval=50.;rmax_use=75000.
    interpmethod='Rbf';data_unrealible=15.0;maxthreshold=3. # maxthreshold = liat dipaper hess-13-195-2009 >> hasil koreksi lenih dari 3 kali dari hasil non koreksi
##    gauge_not_use_opr_blck=read_gauge_not_use(filegaugenotuse)
    sitecoords,size_grid,x_grid,y_grid=sizegrid('{}/{}'.format(radarpath,os.listdir(radarpath)[0]),res,rmax_use)
#    data,sitecoords,xgrid,ygrid,xymaxmin=EEC('{}/{}'.format(radarpath,os.listdir(radarpath)[0]),res,rmax_use,nelevation)
##    gauge_use=read_gauge_use(filemetadata,provinsi,gauge_not_use_opr_blck,sitecoords,rmax_use)

#20190326
##    gauge_use=['ARG Kertajati','AWS Stamet Jatiwangi','ARG Subang','ARG Pusaka Negara','ARG Pabrik Gula Subang','ARG Jampang Kulon','AWS Cisolok',\
##               'AWS Cimalaka','ARG Singaparna','ARG Salopa','AWS Tasikmalaya','AWS Stageof Bandung','AWS IPB','ARG Rekayasa Cisadane',\
##               'ARG Kebun Raya Bogor','AWS UI','ARG Cidaun','ARG Sukanegara','ARG Ciranjang',\
##               'ARG Banjar Irigasi','AAWS Lebak',\
##               'AWS SMPK Cileles','AWS Ujung Kulon','ARG Menes','ARG Cisalak','ARG Cibaliung','AWS Stamet Serang','ARG Ciomas',\
##               'ARG Cikeusal Timur','AWS Pelabuhan Serang','ARG SMPK Singamerta','Stasiun Meteorologi Budiarto Curug','AWS BSD Serpong','AWS Golf Modern','AWS Maritim Merak',\
##               'ARG Kramat Watu','AWS Digi Stamet Serang','ARG Tirtayasa','AWS Digi Stamet Cengkareng','ARWS Rekayasa Stageof Tangerang','ARG BPP Balaraja Tangerang','Stasiun Klimatologi Tangerang Selatan',\
##               'BPP Mauk','AWS Staklim Tangerang Selatan','ARG Rekayasa STMKG',\
##               'AWS Cikancung (Ex AWS Solokan Jeruk)','ARG Bekasi','AAWS Bojongmangu','ARG Cariu',\
##               'ARG Cikasungka','AWS Puspitek','Stasiun Klimatologi Bogor','AAWS Dramaga','AWS Leuwiliang','AWS Mekarsari','AWS Jagorawi',\
##               'AWS Gunung Geulis','AWS Parung','AWS Stamet Citeko','AWS Cibeureum','ARWS Rekayasa Cibeureum',\
##               'ARG Rekayasa Pintu Air Cibongas','ARG Rekayasa Citeko','ARG Rekayasa Bendungan Cipamingkis','ARG Kawali','AAWS Banjarsari Ciamis','ARG Gegesik','ARG Setu Patok',\
##               'ARG Ciberes']

#20181212
    gauge_use=['AWS Stamet 745 Kemayoran','ARG Rekayasa Kemayoran','ARG Manggarai','ARG Lebak Bulus','AWS TMII',\
               'AWS Cikancung (Ex AWS Solokan Jeruk)','ARG Bekasi','AAWS Bojongmangu','ARG Cariu','ARG Cikasungka','AWS Puspitek','Stasiun Klimatologi Bogor',\
               'AAWS Dramaga','AWS Leuwiliang','AWS Mekarsari','AWS Jagorawi','AWS Gunung Geulis','AWS Parung','AWS Stamet Citeko',\
               'AWS Cibeureum','ARWS Rekayasa Cibeureum','ARG Rekayasa Pintu Air Cibongas','ARG Rekayasa Citeko','ARG Rekayasa Bendungan Cipamingkis',\
               'ARG Kawali','AAWS Banjarsari Ciamis','ARG Gegesik','ARG Setu Patok','ARG Ciberes','AAWS Lemah Abang',\
               'ARG Cisompet','AWS Cisurupan','AAWS Indramayu 2 (AWS )','AWS Losarang','AAWS Indramayu','AAWS Kota Baru (AWS )',\
               'ARG Rengasdengklok','ARG Karawang','AWS Kadugede','ARG Sukahaji','ARG Kertajati','AWS Stamet Jatiwangi','AWS Digi Stamet Kertajati',\
               'ARG Subang','ARG Pabrik Gula Subang','ARG Jampang Kulon','AWS Cisolok','AAWS Sumedang','ARG Salopa','AWS STT Cipasung',\
               'AWS Tasikmalaya','ARG Ciwidey','AWS Stageof Bandung','AWS IPB','ARG Rekayasa Cisadane','ARG Kebun Raya Bogor','AWS UI','ARG SMPK Tasikmalaya','ARG Cidaun',\
               'ARG Sukanegara','ARG Agrabinta','ARG Ciranjang',\
               'ARG Malingping','ARG Banjar Irigasi','ARG Bojong Leles','AWS Leuwidamar','AWS SMPK Cileles','ARG Menes','ARG Cisalak','ARG Pandeglang',\
               'ARG Cibaliung','AWS Labuhan ','ARG Padarincang','ARG Ciomas','ARG Cikeusal Timur','AAWS Serang','ARG SMPK Singamerta','AWS BSD Serpong',\
               'AWS Golf Modern','AWS Maritim Merak','ARG Kramat Watu','AWS Digi Stamet Serang','ARG Tirtayasa','AWS Digi Stamet Cengkareng',\
               'ARWS Rekayasa Stageof Tangerang','ARG BPP Balaraja Tangerang','AWS Staklim Tangerang Selatan','ARG Rekayasa STMKG']
               #,'ARG Kelapa Gading']


# lat,lon,qpe=loadnc(fname)
#convert to csv
#df = pd.DataFrame(radar_acc)
#df.to_csv('radar_acc_2019032613.csv', sep=',', index=False)
# interpmethod dapat diganti data dengan : griddata atau Rbf
# pilihan griddata  : nearest; linear; cubic
# pilihan Rbf       : cubic, gaussian, inverse, linear, quintic, thin_plate


    try:
        sfpath,stime=collectdata(radarpath,timeacc,vendor)
        data_qpe=calculate_qpe(radarpath,sfpath,stime,res,nelevation,size_grid,rmax_use)
        gauge_acc,radar_acc,gauge_lon,gauge_lat=datapoint(filemetadata,fileaws,provinsi,gauge_use,y_grid,x_grid,data_qpe,timeacc,data_unrealible)
        print (len(gauge_use))
        print (len(gauge_acc))
        print (len(radar_acc))
        print(radar_acc)
        print (c)
        correction_mfb,correction_bra,correction_add,grid_mfb,grid_bra,grid_add,grid_gauge=adjust(gauge_acc,radar_acc,gauge_lon,gauge_lat,sitecoords,rmax_use,x_grid,y_grid,interpmethod,maxthreshold)
        write_csv(gauge_use,gauge_acc,radar_acc,timeacc,csvpath,correction_mfb,correction_bra,correction_add)
        mae,me,rmse,mae_adj,me_adj,rmse_adj,bias=verification(gauge_acc,radar_acc,correction_mfb,adjtype='multiplicative')
        mae,me,rmse,mae_adj,me_adj,rmse_adj,bias=verification(gauge_acc,radar_acc,correction_bra,adjtype='multiplicative')
        mae,me,rmse,mae_adj,me_adj,rmse_adj,bias=verification(gauge_acc,radar_acc,correction_add,adjtype='additive')
        title='Radar Unadjusted '+timeacc.strftime("(%Y%m%d%H00)");file1='{}/{}_{}.png'.format(imagepath,title,timeacc.strftime("%Y%m%d%H00"))
        plotonmap(y_grid,x_grid,shppath,data_qpe,title,file1,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type='rain')
        title='Mean Field Bias Adjustment '+timeacc.strftime("(%Y%m%d%H00)");file1='{}/{}_{}_{}.png'.format(imagepath,title,timeacc.strftime("%Y%m%d%H00"),interpmethod)
        plotonmap(y_grid,x_grid,shppath,data_qpe*grid_mfb,title,file1,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type='rain')
        title='Brandes Spatial Adjustment '+timeacc.strftime("(%Y%m%d%H00)");file1='{}/{}_{}_{}.png'.format(imagepath,title,timeacc.strftime("%Y%m%d%H00"),interpmethod)
        plotonmap(y_grid,x_grid,shppath,data_qpe*grid_bra,title,file1,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type='rain')
        title='Additive Spatial Adjustment '+timeacc.strftime("(%Y%m%d%H00)");file1='{}/{}_{}_{}.png'.format(imagepath,title,timeacc.strftime("%Y%m%d%H00"),interpmethod)
        plotonmap(y_grid,x_grid,shppath,data_qpe+grid_bra,title,file1,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type='rain')
        title='Gauge Interpolation '+timeacc.strftime("(%Y%m%d%H00)");file1='{}/{}_{}_{}.png'.format(imagepath,title,timeacc.strftime("%Y%m%d%H00"),interpmethod)
        plotonmap(y_grid,x_grid,shppath,grid_gauge,title,file1,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type='rain')
        title='MFB Correction Factor '+timeacc.strftime("(%Y%m%d%H00)");file1='{}/{}_{}_{}.png'.format(imagepath,title,timeacc.strftime("%Y%m%d%H00"),interpmethod)
        plotonmap(y_grid,x_grid,shppath,grid_mfb,title,file1,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type='correction')
        title='BRA Correction Factor '+timeacc.strftime("(%Y%m%d%H00)");file1='{}/{}_{}_{}.png'.format(imagepath,title,timeacc.strftime("%Y%m%d%H00"),interpmethod)
        plotonmap(y_grid,x_grid,shppath,grid_bra,title,file1,gauge_acc,correction_mfb,correction_bra,gauge_lon,gauge_lat,type='correction')
    except:errlog()


#df = pd.DataFrame(radar_acc)
#df.to_csv('radar_acc_2019032613.csv', sep=',', index=False)
