def flatten(lst):
    a = [item for sublist in lst for item in sublist]
    return a

def movmean(lst,n,step,dim):
    # moving average

    import numpy as np

    try:
        if np.ndim(lst)>1:
            if dim == 0:
                a = [[np.nanmean(lst[np.arange(0, n) + i, j]) for i in range(0, np.size(lst, 0) - n + 1, step)]
                     for j in range(0, np.size(lst, 1))]
            elif dim == 1:
                a = [[np.nanmean(lst[j, np.arange(0, n)+i]) for i in range(0, np.size(lst, 1) - n + 1, step)]
                     for j in range(0, np.size(lst, 0))]
        else:
            a = [np.nanmean(lst[np.arange(0, n)+i]) for i in range(0, len(lst)-n+1, step)]
    except:
        n

    return a

def peaksearch(lst,threshold,option):
    # search peak that is above threshold

    # lst: signal
    # threshold: threshold for peak
    # option: the way you find peak - find max, first, or last peak
    import numpy as np
    if len(threshold) == 1:
        threshold = np.repeat(threshold, len(lst))

    peakidx = [[i for i, v in enumerate(lst[j]) if v > threshold[j]] for j in range(0, np.size(lst, 0))]
    peakval = [[v for i, v in enumerate(lst[j]) if v > threshold[j]] for j in range(0, np.size(lst, 0))]

    if option != 'nan':
        if option == 'max':
            peakidx = [x[np.argmax(y)] if len(x) > 0 else np.nan for x, y in zip(peakidx, peakval)]
            peakval = [np.max(x) if len(x)>0 else np.nan for x in peakval]
        elif option == 'first':
            peakidx = [x[0] if len(x) > 0 else np.nan for x in peakidx]
            peakval = [x[0] if len(x) > 0 else np.nan for x in peakval]
        elif option == 'last':
            peakidx = [x[-1] if len(x) > 0 else np.nan for x in peakidx]
            peakval = [x[-1] if len(x) > 0 else np.nan for x in peakval]

    return peakidx, peakval


def lnregress(X, y):
    # linear regression
    from scipy import stats
    from sklearn.linear_model import LinearRegression
    import numpy as np
    lm = LinearRegression()
    lm.fit(X, y)
    params = np.append(lm.intercept_, lm.coef_)
    predictions = lm.predict(X)
    new_X = np.append(np.ones((len(X), 1)), X, axis=1)
    M_S_E = (sum((y - predictions) ** 2)) / (len(new_X) - len(new_X[0]))
    v_b = M_S_E * (np.linalg.inv(np.dot(new_X.T, new_X)).diagonal())
    s_b = np.sqrt(v_b)
    tstat = params / s_b
    p_val = [2 * (1 - stats.t.cdf(np.abs(i), (len(new_X) - len(new_X[0])))) for i in tstat]

    return p_val, tstat, params