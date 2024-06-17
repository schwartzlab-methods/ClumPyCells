import math
import copy
import sys,logging
import numpy as np


def closepairs(X, rmax, d=None, twice=True, what=None,
                           distinct=True, neat=True, periodic=True, pp=None):
    if not rmax >= 0:
        sys.exit(2)
    if not X.getWindow().getType()=="rectangular":
        logging.warning("Periodic edge correction applied in non-rectangular window")
    npts = len(X.getX())
    if not what:
        what = "all"
    if what == "all":
        null_answer = {"i":0,"j":0,"xi":0,
                    "yi":0,"xj":0,"yj":0,
                    "dx":0,"dy":0,"d":0}
    elif what == "indices":
        null_answer = {"i":0,"j":0}
    elif what == "ijd":
        null_answer = {"i":0,"j":0,"d":0}
    if npts==0:
        return null_answer
    
    xsort = copy.deepcopy(X.getX())
    xsort.sort()
    ysort = copy.deepcopy(X.getY())
    ysort.sort()
    # calculate a conservative estimate of the number of pairs
    npairs = math.pow(float(npts),2)
    if(npairs <= 1024):
        nsize = 1024
    else:
        catchfraction = math.pi * math.pow(rmax,2)/X.getWindow().getArea()
        nsize = math.ceil(2 * catchfraction * npairs)
        nsize = min(nsize, npairs)
        nsize = max(1024, nsize)
        
    # Now extract pairs
    if periodic:
        # special algorithm for periodic distance
        got_twice = True
        x = X.getX()
        y = X.getY()
        d = X.getD()
        if d is None:
            d = [0]*len(x)
        r = rmax
        p = X.getSideLength()
        ng = nsize
        z = closePpairs(xx=x,yy=y,dia=d,rr=r,nguess=ng,pp=pp)
        i = z[0]
        j = z[1]
        d = z[2]
        areaWt = z[3]
        if (what == "all"):
            xi = [x[index-1] for index in i]
            yi = [y[index-1] for index in i]
            xj = [x[index-1] for index in j]
            yj = [y[index-1] for index in j]
            dx = [element1 - element2 for (element1, element2) in zip(xj, xi)]
            dy = [element1 - element2 for (element1, element2) in zip(yj, yi)]
    else:
        if not distinct:
            ii = list(range(npts))
            xx = X.getX()
            yy = X.getY()
            zeroes = [0]*npts
            if what == "all":
                null_answer = {"i":ii,"j":ii,"xi":xx,"yi":yy,
                                "xj":xx,"yj":yy,"dx":zeroes,
                                "dy":zeroes,"d":zeroes}
            elif what == "indices":
                null_answer = {"i":ii,"j":ii}
            elif what == "ijd":
                null_answer = {"i":ii,"j":ii,"d":zeroes}
        got_twice = True
        z = Fclosepairs(npts,xsort,ysort,rmaxplus)
        if z["status"] != 0:
            rmaxplus = 1.25 * rmax
            nsize = paircount(npts,xsort,ysort,rmaxplus)
            if nsize<=0:
                return null_answer
            # add a bit more for safety
            nsize = math.ceil(1.1 * nsize) + 2 * npts
            # now extract points
            z = Fclosepairs(npts,xsort,ysort,rmaxplus)

        # trim vectors to the length indicated
        npairs = z["nout"]
        if(npairs <= 0):
            return null_answer
        actual = list(range(npairs))
        i = [z[0][index] for index in actual]
        j = [z[1][index] for index in actual]
        if what == "all":
            xi = [z[2][index] for index in actual]
            xj = [z[3][index] for index in actual]
            yi = [z[4][index] for index in actual]
            yj = [z[5][index] for index in actual]
            dx = [z[6][index] for index in actual]
            dy = [z[7][index] for index in actual]
            d = [z[8][index] for index in actual]
        elif what == "ijd":
            d = [z[8][index] for index in actual]

    if twice:
        if not got_twice:
            # duplication required
            iold = i
            jold = j
            i = iold+jold
            j = jold+iold
            if what == "all":
                xinew = xi+xj
                yinew = yi+yj
                xjnew = xj+xi
                yjnew = yj+yi
                xi = xinew
                yi = yinew
                xj = xjnew
                yj = yjnew
                dx = dx+[ -x for x in dx]
                dy = dy+[ -y for y in dy]
                d = [d, d]

            elif what == "ijd":
                d = [d,d]
    else:
        if got_twice:
            # remove duplication
            if (i<j):
                ok=1
            else:
                ok=0
            i=i[ok]
            j=j[ok]
            if what == "all":
                xi = xi[ok]
                yi = yi[ok]
                xj = xj[ok]
                yj = yj[ok]
                dx = dx[ok]
                dy = dy[ok]
                d  =  d[ok]
            elif what == "ijd":
                d  =  d[ok]
        elif neat:
            # enforce i < j
            if (i>j):
                swap=1
            else:
                swap=0
            tmp = i[swap]
            i[swap] = j[swap]
            j[swap] = tmp
            if what == "all":
                if swap==1:
                    xinew = xj
                    xjnew = xi
                    yinew = yj
                    yjnew = yi
                else:
                    xinew = xi
                    xjnew = xj
                    yinew = yi
                    yjnew = yj
            xi = xinew
            yi = yinew
            xj = xjnew
            yj = yjnew
            dx[swap] = -dx[swap]
            dy[swap] = -dy[swap]
    if not distinct:
        ii = list(range(npts))
        xx = X.getX()
        yy = X.getY()
        zeroes = [0]*npts
        i = i+ii
        j = j+ii
        if what =="all":
            xi = xi+xx
            yi = yi+yy
            xj = xj+xx
            yj = yj+yy
            dx = dx+zeroes
            dy = dy+zeroes
            d  = d+zeroes
        elif what =="ijd":
            d = d+zeroes
    if what=="all":
        answer = {"i":i,"j":j,"xi":xi,
                    "yi":yi,"xj":xj,"yj":yj,
                    "dx":dx,"dy":dy,"d":d, "Awt" : areaWt}
    elif what=="indices":
        answer = {"i":i,"j":j}
    elif what=="ijd":
        answer = {"i":i,"j":j,"d":d}
    return answer


def paircount(nxy, x, y, rmaxi):
    # nxy: number of (x,y) points
    # x, y: (x,y) coordinates
    # rmaxi: maximum distance
    count = 0
    r2max = rmaxi * rmaxi
    if nxy==0:
        return
    i = 0
    maxchunk = 0
    while (i<nxy):
        maxchunk += 65536
        if(maxchunk > nxy):
            maxchunk =nxy
        while i<maxchunk:
            xi = x[i]
            yi = y[i]
            if (i>0):
                for j in range(i - 1,-1,-1):
                    dx = x[j] - xi
                    a = r2max - dx * dx
                    if(a < 0):
                        break
                    dy = y[j] - yi
                    a -= dy * dy
                    if(a >= 0):
                        count=count+1
            if (i+1<nxy):
                for j in range(i+1, nxy, 1):
                    dx = x[j] - xi
                    a = r2max - dx * dx
                    if(a < 0):
                        break
                    dy = y[j] - yi
                    a -= dy * dy
                    if(a >= 0):
                        count=count+1
            i=i+1
    return count


def closePpairs(xx, yy, dia, rr, nguess, pp=None):
    # pp is the point pattern set for adipocyte
    n = len(xx)
    r2max = rr*rr
    pp_x = []
    pp_y = []
    pp_d = []
    if pp is not None:
        pp_x = pp.getX()
        pp_y = pp.getY()
        pp_d = pp.getD()
    jout, iout, dout=[],[],[]
    areaWt = []
    if(n > 0 and nguess > 0):
        i = 0
        maxchunk = 0
        while(i < n):
            maxchunk += 65536
            if (maxchunk > n):
                maxchunk = n
            while i<maxchunk:
                xi = xx[i]
                yi = yy[i]
                di = dia[i]
                if (i>0):
                    for j in range(i - 1,-1,-1):
                        dx = xx[j] - xi
                        if(dx < 0.0):
                            dx = -dx
                        if(dx < rr):
                            dy = yy[j] - yi
                            if(dy < 0.0):
                                dy = -dy                       
                            d2 = max(0, math.sqrt(dx * dx + dy * dy) - di/2 - dia[j]/2)
                            distance = math.sqrt(d2)
                            area = overlapA(di, distance, distance)
                            for k in range(len(pp_x)):
                                p1 = np.array([xi,yi])
                                p2 = np.array([xx[j],yy[j]])
                                p3 = np.array([pp_x[k],pp_y[k]])
                                d = np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1) # Distance between point and a line (from two points)
                                if d<pp_d[k]/2:
                                    dist_pass = 2*math.sqrt(pp_d[k]**2/4-d**2)
                                    if dist_pass < d2:
                                        d2 = d2 - dist_pass
                            if(d2 * d2 <= r2max):
                                jout.append(j + 1)
                                iout.append(i + 1)
                                dout.append(d2)
                                areaWt.append(area)
                if(i + 1 < n):
                    for j in range(i+1,n,1):
                        dx = xx[j] - xi
                        if(dx < 0.0):
                            dx = -dx
                        if(dx < rr):
                            dy = yy[j] - yi
                            if(dy < 0.0):
                                dy = -dy
                            d2 = max(0, math.sqrt(dx * dx + dy * dy) - di/2 - dia[j]/2)
                            distance = math.sqrt(d2)
                            area = overlapA(di, distance, distance)
                            for k in range(len(pp_x)):
                                p1 = np.array([xi,yi])
                                p2 = np.array([xx[j],yy[j]])
                                p3 = np.array([pp_x[k],pp_y[k]])
                                d = np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1) # Distance between point and a line (from two points)
                                if d<pp_d[k]/2:
                                    dist_pass = 2*math.sqrt(pp_d[k]**2/4-d**2)
                                    if dist_pass < d2:
                                        d2 = d2 - dist_pass
                                    # print("d2=",d2)
                            if(d2 * d2 <= r2max):
                                jout.append(j + 1)
                                iout.append(i + 1)
                                dout.append(d2)
                                areaWt.append(area)
                i = i+1
    return iout,jout,dout,areaWt


def overlapA (r1, r2, d):
    # This function returns the overlap area of two cicles, r1 is the size of the point (cell), r2 is the radius of the disk, d is the distance between the two points
    if d >= r1 + r2:
        return 0
    
    elif d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2

    alpha1 = math.acos((d ** 2 + r1 ** 2 - r2 ** 2) / (2 * d * r1))
    alpha2 = math.acos((d ** 2 + r2 ** 2 - r1 ** 2) / (2 * d * r2))
    
    A1 = r1 ** 2 * alpha1
    A2 = r2 ** 2 * alpha2
    A3 = -0.5 * ((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)) ** 0.5
    
    return A1 + A2 + A3


def Fclosepairs(nxy, x, y, r, noutmax):
    r2max = r*r
    jout,iout,xjout,xiout,yjout,yiout,dxout,dyout,dout=[]
    if (nxy==0):
        return
    i = 0
    k = 0
    maxchunk = 0
    while (i<nxy):
        maxchunk += 65536
        if (maxchunk > nxy):
            maxchunk = nxy
        while i<maxchunk:
            xi = x[i]
            yi = y[i]
            if (i>0):
                for j in range(i - 1,-1,-1):
                    dx = x[j] - xi
                    dx2 = dx * dx
                    if(dx2 > r2max):
                        break
                    dy = y[j] - yi
                    d2 = dx2 + dy * dy
                    if(d2 <= r2max):
                        if (k>=noutmax):
                            return
                        jout.append(j + 1)
                        iout.append(i + 1)
                        xiout.append(xi)
                        yiout.append(yi)
                        xjout.append(x[j])
                        yjout.append(y[j])
                        dxout.append(dx)
                        dyout.append(dy)
                        dout.append(math.sqrt(d2))
                        k=k+1
            if (i+1<nxy):
                for j in range(i+1, nxy, 1):
                    dx = x[j] - xi
                    dx2 = dx * dx
                    if(dx2 > r2max):
                        break
                    dy = y[j] - yi
                    d2 = dx2 + dy * dy
                    if(d2 <= r2max):
                        if (k>=noutmax):
                            return
                        jout.append(j + 1)
                        iout.append(i + 1)
                        xiout.append(xi)
                        yiout.append(yi)
                        xjout.append(x[j])
                        yjout.append(y[j])
                        dxout.append(dx)
                        dyout.append(dy)
                        dout.append(math.sqrt(d2))
                        k=k+1
            i=i+1
    return iout,jout,xiout,xjout,yiout,yjout,dxout,dyout,dout,k
