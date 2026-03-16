import pandas as pd
import numpy as np




# Helper function to perform logarithmic interpolation of missing points between minimum values
def hysep_interp(df):
    nf = 10 ** np.log10(df).interpolate(method='linear')
    return nf


# Baseflow separation program based on USGS HYSEP program
# local minimum method only
# USGS Report for Hydrograph Separation routine can be found at:
# http://www.water-research.net/Waterlibrary/geologicdata/estbaseflow.pdf
#
# datain = n x 2
#    datain(:,1) = timestep (daily)
#    datain(:,2) = flow rate (mean daily)
# Area = square miles
#
# Output:
#  dataout = n x 4 [time baseflow stormflow totalflow]
#
# Syntax:
#   dataout = f_hysep(datain, Area);
#
def hysep(df, area):
    # The duration of surface runoff is calculated from the empirical relation:
    # N=A^0.2, (1) where N is the number of days after which surface runoff ceases,
    # and A is the drainage area in square miles (Linsley and others, 1982, p. 210).
    # The interval 2N* used for hydrograph separations is the odd integer between
    # 3 and 11 nearest to 2N (Pettyjohn and Henning, 1979, p. 31).
    N = area ** 0.2
    
    # The hydrograph separation begins one interval (2N* days) prior to the start of the date selected
    # for the start of the separation and ends one interval (2N* days) after the end of the selected date
    # to improve accuracy at the beginning and end of the separation. If the selected beginning and
    # (or) ending date coincides with the start and (or) end of the period of record, then the start of the
    # separation coincides with the start of the period of record, and (or) the end of the separation
    # coincides with the end of the period of record.
    if np.ceil(2*N) % 2 == 0:
        inN = np.floor(2*N)
    else:
        inN = np.ceil(2*N)
    inN = int(min(max(inN, 3), 11))

        
    #cf = df.rolling(window=wdw, center=True).min()
    #cf.loc[(df.Flow > cf.Flow),'Flow'] = np.nan

    # Filter out lowest flow values per a rolling window
    cf = df.rolling(window=inN, center=True, min_periods=0).agg(np.nanmin)
    cf['Flow'] = (df.Flow <= cf.Flow)*df.Flow
    cf = cf.replace(0, np.nan)
    cf[:int(0.5*(inN-1))] = np.nan
    cf[-int(0.5*(inN-1)):] = np.nan
    first_valid = cf[cf.Flow.notnull()].index[0]
    last_valid = cf[cf.Flow.notnull()].index[-1]
    cf = cf[first_valid:last_valid]
    

    # Interpolate the missing values and replace values higher than original
    cf['minFlow'] = hysep_interp(cf[['Flow']])
    cf['Flow'] = cf.min(axis=1)
    cf = cf[['Flow']]

    cf.join(df, lsuffix='total', rsuffix='base')
    mf = cf.join(df, how='right', lsuffix='base', rsuffix='total')
    mf['Flowstream'] = round(mf.Flowtotal - mf.Flowbase, 2)
    mf = mf[[first_valid:last_valid]]
    mf = mf[['Flowbase', 'Flowstream', 'Flowtotal']]
    return mf




# stormevents - find and extract storm events in a timeseries
# This function will identify and extract storm event flow rates using HYDSEP
# techniques (local minimum).  The user needs to provide a time series of data,
# the drainage area associated with the time series, a set of peak flow rates
# to use, and a tolerance of the peak flow rate to be used. The original
# usage of HYDSEP was meant to be used on daily mean flow rates, this
# function uses maximum daily flow rates instead.
def stormevents(df, a, q_lvl, q_tol):
    # get daily max
    df = df.resample('1d').max()

    hysep1 = hysep(df, a)

    s = (hysep1.Flowstream() > 0)*1
    

    maxevt = [hysep1()]
    s = hysep1()
    u = (s - s.shift(1)).dropna()
    u = u.replace(-1,0)
    v = u.cumsum()
    hysep1['evt'] = v
    

    
    
    
    s = hysep1(:,3) > 0;
    t = double(s);
    u = diff(t);
    u(u == -1) = 0; #remove negative flag for ending an event, don't need.
    v = cumsum(u); %event number
    hysep1 = [v hysep1(1:end-1,:)]; % insert event number as first column

    clear s t u

    maxevt = accumTS([v hysep1(:,end)],@max); %roll up peak q to the event number
    [~, ia, ~] = intersect([hysep1(:,1) hysep1(:,5)],maxevt,'rows'); % locate dates of event peaks
    maxevt = [hysep1(ia,2) maxevt]; %day of peak q per event

    evts = hysep1(hysep1(:,4) > 0,:); % remove non-event days

    clear hysep1

    qd = q_levels * q_tolerance;
    qs = repmat(maxevt(:,3),1,length(q_levels));
    ql = repmat(q_levels,length(maxevt(:,3)),1);
    qdm = repmat(qd,length(maxevt(:,3)),1);
    d = abs(ql - qs);
    idx = d <= qdm;

    clear qd qs ql qdm d

    [y2, m2, d2, h2, mn2, s2] = datevec(c(:,1));
    idxx = zeros(length(y2),length(q_levels));
    for iq = 1:length(q_levels)
        qevts =  maxevt(idx(:,iq),:);
        iqm = ismember(evts(:,1),qevts);
        [y, m, d] = datevec(evts(iqm,2));
        idx2 = ismember([y2 m2 d2],[y m d],'rows');
        idxx(:,iq) = idx2;
        q{iq} = [y2(idx2) m2(idx2) d2(idx2) h2(idx2) mn2(idx2) s2(idx2) c(idx2,2)];
    end

end
