#import qt
import numpy as np
import scipy.constants as co
from scipy import signal
import random as rnd

from matplotlib import pyplot as plt

gamma = -2*np.pi*28.0e9
BperI = .0698
#TODO use a single FperI parameter instead
d = 8.5e-6
Ms = 197.4e3 #at 10mK
wM = - gamma * co.mu_0 * Ms
ad = .0060 #interantenna distance
pi = np.pi
aw = 50e-6 #antenna width

def freq_to_current(freq):
    '''
    get the current needed to go to an FMR frequency
    Input: frequency in GHz
    Output: current (A)
    '''
    #freq = offset + current*bperI*gamma

    #This is based on the dispersion relation below. Relations inverted using mathematica file _users/arjan/mathematica
    
    w = 2 * np.pi * freq
    I = (wM + np.sqrt(4*w*w + wM*wM)) / (2*BperI*gamma)
    return I

def current_to_freq(I):
    ''' 
    returns the FMR mode for a magnetic field due to a current I

    Input: current (A)
    Output: frequency (Hz)
    '''
    w = np.sqrt(wHI(I)*(wHI(I) + wM))
    return w / (2*np.pi)


##########################################################################################
#---------------------FUNCTIONS FOR SURFACE WAVES----------------------------------------
##########################################################################################


def w_k(I, k):
    '''
    dispersion for surface travelling waves
    '''
    return np.sqrt(wHI(I)*(wHI(I)+wM) + (wM*wM/4)*(1 - np.exp(-2*k*d)))

def w_f(I, f):
    return w_k(I, k_f(I, f))

def wHI(I):
    '''
    omega H as a function of current
    '''
    return -1 * gamma * BperI * I

def f0(I):
    '''
    Get the lowest frequency, associated with the FMR
    '''

    return w_k(I, 0) / (2*np.pi)


def wH(B):
    return -1 * gamma * B

def vg_k(I, k):
    '''
    This is the non-approximated version of dw/dk
    '''
    top = d * np.exp(-2*d*k) * wM
    bottom = 4 * np.sqrt((1./4) * 1 - np.exp(-2*d*k) * wM**2 + wHI(I) * (wHI(I) + wM))
    
    return top / bottom

def k_f(I, f):
    w = 2*np.pi*f
    top = (2*wHI(I) - 2*w + wM)*( 2* (wHI(I) + w) + wM)
    if top > 0:
        k = -1./(2*d) * np.log( top / (wM*wM) )
        return k
    else: 
        return -1


def vg_f(I, f):
    w = 2 * np.pi * f
    vg = d * ((2*wHI(I)+ wM)**2 - 4*w**2) / (4*np.abs(w))
    return vg

def ta_f(I, f):
    '''
    returns the arrival time if the frequency is above f0, otherwise it returns 0
    '''
    print 'static f is ', f
    print f, f0(I)
    return (ad / vg_f(I, f)) * ( f >= f0(I))


def vPhase(I, f):
    if k_f(I,f) <= 0:
        return 0
    else:
        return (vg_f(I, f)/ k_f(I,f))


def pulseDelay(I, f, d=ad):
    return d/vg_f(I,f)


def lorentzian(x, G):
    return (1./pi) * (G/2) / (x*x + (G*G/4))


def RLCresponse(times, R, L, C, B1, B2, plot=False):
    '''
    SERIES RLC
    R = .5
    L = 1e-9
    C = 1e-9
    B1 = .2
    B2 = 1
    '''
    R = float(R)
    alpha = R/(2*L)
    damp = R/2. * np.sqrt(C/L) #damping factor: underdamped if smaller than 1
    print 'damping factor is ', damp
    wr = 1./np.sqrt(L*C) #resonance
    print 'resonance freq is ', wr/(2*pi)/1e9, 'GHz'
    wd = wr * np.sqrt(1-damp*damp)
    
    if damp < 1:
        def resp(t): 
            return B1*np.exp(-alpha*t)*np.cos(wd*t) + B2 * np.exp(-alpha*t) * np.sin(wd*t)
    elif damp > 1:
        print 'overdamped regime does not work yet!!'
        def resp(t):
            return .2*np.exp(damp + np.sqrt(damp*damp - 1)*-1*wr*t) + .2*np.exp(damp - np.sqrt(damp*damp - 1)*-1*wr*t) 
    elif damp == 1:
        print 'critically damped!'
        def resp(t):
            return B1 * t * np.exp(-alpha*t) + B2 * np.exp(-alpha*t)

    
    if plot:
        plotr2(times, resp(times))
        q = 1./R * np.sqrt(L/C)
       # print 'Q is ', Q
    return resp(times)



def CopperStripInd():
    '''
    http://www.k7mem.com/Electronic_Notebook/inductors/wire_strip.html
    '''
    w=50e-6 #width
    b=10e-3 #length
    K=16e-6 #thickness of copper???? This is a guess
    h=.5e-3 #Thickness of material

    L = K*b*(np.log(2*b/(w+h)) + .5 + .2235*(w+h)/b)
        

#=====================================================================================================================
#---------------------------------------------BLOB SIMULATION---------------------------------------------------------
#=====================================================================================================================

pulseStart = 25e-9
pulseLen = 25e-9
pulseAmp = 1
f0 = 7.2e9
w0 = 2*np.pi * f0
IF = 250e6
LO = f0 - IF

def squarePulse(t, dc = False, offset=0):
    if type(t) == np.ndarray:
        k = [int(pulseStart < ti < (pulseStart + pulseLen)) + offset for ti in t]
        return k * pulseAmp
    elif (type(t) == int) or (type(t) == float):
        return int(pulseStart < t < (pulseStart+pulseLen))

def squarePulseSin(t, offset=0):
    if type(t) == np.ndarray:
        k = [np.sin(w0*ti) * (int(pulseStart < ti < (pulseStart + pulseLen)) + offset) for ti in t]
        return k * pulseAmp
    elif (type(t) == int) or (type(t) == float):
        return int(pulseStart < t < (pulseStart+pulseLen))

def gaussPulse(t, DC=False):
    if DC==False:
        k = np.sin(w0*t)
    else:
        k = np.ones(len(t))
    pulseCent = pulseStart + pulseLen/2
    width = pulseLen/4.3
    k = k * np.exp(-(t-pulseCent)**2/(2*width**2))
    plotr2(t, k)
    return k*pulseAmp

def couplingk(k):
    A = 1.2 #random coupling enhancer
    return (A * .5*(1 + np.cos(2*k*aw))) * int( k < .5* pi/aw)

def coupling(I, f):
    k = k_f(I,f)
    if k > 0:
        return couplingk(k)
    return 0



#Mix the input pulse with the antenna BW, which we assume to be lorentzian. We do this by inverse fourier transforming the lorentzian, and convolving the result with our input
def inputAnt(tx, amps, BW=400e6, type='lor', plot=False):
    #time and frequency axis
    dt = tx[1]-tx[0]
    freqs = np.arange(-.5/dt, .5/dt+.001, 1./(tx[-1]-tx[0]))
    
    if type == 'lor':
        #lorentzian peakshape
        antPeak = lorentzian(freqs, BW)
        #calc the inverse FFT of the lorentzian, and shift it to the center
        lortime = np.abs(np.fft.ifftshift(np.fft.ifft(antPeak)))
        #lortime = np.abs(np.fft.ifft(antPeak))
        #plotr2(tx,lortime)
        antSig = np.convolve(amps, lortime, mode='same')
        antSig = antSig / max(antSig)
    if type == 'rlc':
        print 'rlc type!'
        #antenna parameters
        #R = .6
        #L = .002e-9
        #C = .002e-8
        #B1 = .4
        #B2 = 1
        R = .02
        L = .001e-9
        C = .1e-9
        B1 = .4
        B2 = .1
        
        rlcTime = RLCresponse(tx[:len(tx)/2], R, L, C, B1, B2)
        #displace RLCtime to the center
        rlcTime = np.array( list(np.zeros(len(tx)/2)) + list(rlcTime))
        #plotr(rlcTime)
        antSig = np.convolve(amps, rlcTime, mode='same')
    if plot:
        plotr2(tx, antSig)
    return np.array(antSig)

def sampleBoxResponse(tx, amps, plot=False):
    #parameters
    R = .1
    L = .01e-9
    C = .01e-9
    B1 = 0
    B2 = 1

    #response
    rlcTime = RLCresponse(tx[:len(tx)/2], R, L, C, B1, B2)
    #displace RLCtime to the center
    rlcTime = np.array( list(np.zeros(len(tx)/2)) + list(rlcTime))
    #plotr(rlcTime)
    boxSig = np.convolve(amps, rlcTime, mode='same')
    if plot:
        plotr2(tx, boxSig)


    return np.array(boxSig)

def inYigOld(ti, amps, I, f):
    '''
    Needs a lot more thought
    '''

    #delay time
    td = pulseDelay(I, f)

    #assume infinite bandwidth in the yig
    yigResponse = np.zeros(len(ti)) 
    index = td/(ti[1]-ti[0])
    if (td >0) and (index < len(yigResponse)):
        yigResponse[index]=1
    elif td<0:
        print 'I: ',I
    #plotr2(ti, yigResponse) 
    #set 0 times in the center
    yigResponse = list(np.zeros(len(yigResponse))) + list(yigResponse)

    resp = np.convolve(amps, yigResponse, mode='full')[len(ti):2*len(ti)]

    #plotr2(ti, resp)
    return resp

def inYig3Old(ti, amps, I, f0, d=ad):
    '''
    Tries to find response of YIG to full frequency spectrum of the input signal at full time resolution.
    '''
    #find the frequency spectrum of the input signal
    #plotr(amps)
    inputF = np.fft.fft(amps)
    inputF = np.fft.fftshift(inputF)
    
    #find the corresponding frequencies
    sf = 1./(ti[1]-ti[0])
    freqs = np.linspace(-sf/2.,  sf/2.,len(ti))
    if len(freqs) != len(inputF):
        raise ValueError('Lengths of FFT and frequency vector do not match!: ', len(inputF), len(freqs))

    #determine the YIG passband (in terms of our downconverted frequencies)
    wH = wHI(I)
    #print LO/1e6
    wMin = np.sqrt(wH*(wH+wM))
    wMin -= 2*pi*1000e6
    wMax = wH + 0.5*wM

    output = np.zeros(len(amps))
    newWave = []

    #find the time delay due to the magnon propagation
    for i in range(len(freqs)):
        w = 2*np.pi*freqs[i]
        #check the frequency is in the range of allowed MSSW frequencies
        if ((w > wMin) and (w < wMax)) or ((w>-wMax) and (w<-wMin)):
            #if within the YIG band, apply the appropriate phase shift
            td = pulseDelay(I, freqs[i], d)
            #td = pulseDelay(I, f0, d) #TODO The above line caused the strange behaviour
            newWave.append((inputF[i] * np.exp(-1j* w *td)) ) #* (1-(w-wMin)/1e10))
            #TODO WHAT IS THIS EXTRA FACTOR IN THE PREVIOUS LINE?
        else:
            newWave.append(0)
    newWave = np.array(newWave)
    newWave = np.fft.ifftshift(newWave)
    newWave = np.fft.ifft(newWave)
    #plt.plot(ti, amps)
    #plt.plot(ti, newWave)
    #plt.show()
    
    return newWave

   
    
    
def main(pulse='dc'):

    couplingConstant = .1
    
    times = np.arange(0, 796e-9, .4e-9)
    sf = 1./(times[1]-times[0])
    freqs = np.linspace(-sf/2,sf/2.,len(times))
    
    #generate input pulse
    if pulse == 'squ':
        inp = squarePulseSin(times, noise)
    elif pulse == 'dc':
        inp = squarePulse(times, offset=.02)
    elif pulse == 'gauss':
        inp = gaussPulse(times)
    elif pulse == 'gaussDC':
        inp = gaussPulse(times, DC=True)

    #mix with input antenna
    atype = 'rlc'  #antenna type, can be rlc or log
    inpAnt = inputAnt(times, inp, BW=300e6, type=atype) #BW is ignored for rlc filters
    
    #For direct coupling: Mix with box response and output antenna
    boxSig = sampleBoxResponse(times, inpAnt)
    outpAntD = inputAnt(times, boxSig, type=atype)

    plotr2(times, outpAntD, title='outpAntD')

    #couple into yig
    I0 = 2.2
    currents = np.arange(I0-.1, I0+.2, .001)
    data = []
    for I in currents:
        #input coupling due to antenna shape
        inpFFT = np.fft.fft(inpAnt)
        inpFFT = np.fft.fftshift(inpFFT)
        for i in range(len(freqs)):
            inpFFT[i] = coupling(I, freqs[i]+f0) * inpFFT[i]
        
        #yig response:
        fresp =  inYig(freqs, inpFFT, I, f0, ad)
        outpFFT = fresp

        #output coupling, with antenna shape
        for i in range(len(freqs)):
            outpFFT[i] = outpFFT[i] * coupling(I, freqs[i]+f0)
        
        #back to time domain
        outpFFT = np.fft.ifftshift(outpFFT)
        outpYig = np.fft.ifft(outpFFT)

        #total output: Yig plus direct coupling
        output = outpYig + outpAntD
        
        #testPlots
        if (I >= 2.1499 and I <=2.1509) or (I >= 2.2099 and I <= 2.2109):
            plt.plot(times, inp)
            plt.plot(times, np.fft.ifft(fresp))
            plt.show()
            plt.plot(freqs, np.fft.fftshift(np.fft.fft(inp)))
            plt.plot(freqs, fresp)
            #plt.plot(times, output)
            plt.show()
            #return
        data.append(list(output))

    #plot a subset

    data = np.array(data)
    pdata = data[:,0:800]
    cmap = plt.cm.gnuplot

    rdat = np.real(pdata)
    idat = np.imag(pdata)
    sqdat = np.sqrt(rdat*np.conj(rdat) + idat*np.conj(idat))
    sqmdat = np.sqrt(abs(rdat*np.conj(rdat) - idat*np.conj(idat)))

    fig, ((ax1, ax2),(ax3, ax4), (ax5, ax6)) = plt.subplots(3,2)
    im1 = ax1.imshow(np.transpose(rdat), cmap=cmap, aspect=500000, origin='lower')
    im2 = ax2.imshow(np.transpose(idat), cmap=cmap, aspect=500000, origin='lower')
    im3 = ax3.imshow(np.transpose(sqdat), cmap=cmap, aspect=500000, origin='lower')
    im4 = ax4.imshow(np.transpose(sqmdat), cmap=cmap, aspect=500000, origin='lower')
    im5 = ax5.imshow(np.transpose(abs(rdat)), cmap=cmap, aspect=500000, origin='lower')
    im6 = ax6.imshow(np.transpose(abs(idat)), cmap=cmap, aspect=500000, origin='lower')
    im1.set_extent([currents[0], currents[-1], times[0], times[800]])
    im2.set_extent([currents[0], currents[-1], times[0], times[800]])
    im3.set_extent([currents[0], currents[-1], times[0], times[800]])
    im4.set_extent([currents[0], currents[-1], times[0], times[800]])
    im5.set_extent([currents[0], currents[-1], times[0], times[800]])
    im6.set_extent([currents[0], currents[-1], times[0], times[800]])
    ax1.set_title('real')
    ax2.set_title('imag')
    ax3.set_title('i^2 + q^2')
    ax4.set_title('i^2 - q^2')
    ax5.set_title('abs(real)')
    ax6.set_title('abs(imag)')
    fig.suptitle('shifted version')
    plt.show()

    #return pdata
     
def inYig(freqs, inputF, I, f0, d=ad):
    '''
    Tries to find response of YIG to full frequency spectrum of the input signal at full time resolution.
    '''
    #find the frequency spectrum of the input signal
    sf = freqs[-1] - freqs[0]
    tTime = abs(1./(freqs[1]-freqs[0]))

    #determine the YIG passband (in terms of our downconverted frequencies)
    wH = wHI(I)
    wMin = np.sqrt(wH*(wH+wM))
    #wMin -= 2*pi*1000e6
    wMax = abs(wH + 0.5*wM)

    output = np.zeros(len(inputF))
    newWave = []

    #find the time delay due to the magnon propagation
    for i in range(len(freqs)):
        f = freqs[i]
        w = 2*np.pi*(f + f0)
        #check the frequency is in the range of allowed MSSW frequencies
        xF = 2*pi*1e9
        if (((w >= wMin - xF) and (w < wMax + xF)) or ((w>-wMax-xF) and (w<= -wMin +xF))):
        #if ((w > wMin) and (w < wMax)) or ((w>-wMax) and (w<-wMin)):
            #td = pulseDelay(I, f)#/tTime
            td = pulseDelay(I, f + f0) #TODO The above line caused the strange behaviour
            #TODO Warning! This is cast to real numbers!
            inputF[i] = np.complex(inputF[i]) * np.exp(-1j * td * (f+f0))
            #inputF[i] = inputF[i] * np.real(np.exp(1j * vPhase(I, freqs[i] + f0) * td))
            newWave.append(inputF[i]) #* np.real(np.exp(-1j - vPhase(I, freqs[i]))))#  * (1-(w-wMin)/1e10))
        else:
            newWave.append(0)
    newWave = np.array(newWave)
    
    return newWave

def main3(pulse='squ'):
    '''
    Try to do the simulation more fully - simulate the YIG interaction with much higher time resolution
    and then apply all the measurement-like filtering that we actually use.
    '''

    #couplingConstant = .1
    
    times = np.arange(0, 800e-9, .05e-9)
    sf = 1./(times[1]-times[0])
    print 'sampling frequency', sf/1e9, ' GHz'
    freqs = np.linspace(-sf/2,sf/2.,len(times))
    #generate input pulse
    if pulse == 'squ':
        inp = squarePulseSin(times, offset=.03)
    elif pulse == 'dc':
        inp = squarePulse(times)
    elif pulse == 'gauss':
        inp = gaussPulse(times)

    #mix with input antenna
    atype = 'rlc'  #antenna type, can be rlc or log
    inpAnt = inputAnt(times, inp, BW=300e6, type=atype) #BW is ignored for rlc filters
    #plt.plot(inpAnt)
    #plt.show()
    #plt.plot(freqs, np.fft.rfft(inpAnt))
    #plt.show()
    
    #For direct coupling: Mix with box response and output antenna
    boxSig = sampleBoxResponse(times, inpAnt)
    outpAntD = inputAnt(times, boxSig, type=atype)

    plotr2(times, outpAntD, title='outpAntD')
    #plot fft
    outDF = np.fft.fft(outpAntD)
    outDF = np.fft.fftshift(outDF)
    #print len(outDF)
    tf = np.arange(1, .5/(times[1]-times[0]) + 1, 1./times[-1])
    #print len(tf)
    #plotr(outDF)

    #couple into yig
    I0 = 2.2
    currents = np.arange(I0-.05, I0+.2, .0005)
    data = []
    datar = []
    dataf = []
    
    #examine yig filter
    wH = wHI(I0)
    wMin = np.sqrt(wH*(wH+wM))
    #wMin -= 2*pi*1000e6
    wMax = abs(wH + 0.5*wM)
    print 'wMin is ', wMin/(1e9 * 2 * pi), ' GHz'
    print 'wMax is ', wMax/(1e9 * 2 * pi), ' GHz'

    for I in currents:
        #input coupling due to antenna shape
        inpFFT = np.fft.fft(inpAnt)
        inpFFT = np.fft.fftshift(inpFFT)
        #inpFFT0 = inpFFT.copy()
        for i in range(len(freqs)):
            inpFFT[i] = inpFFT[i] * coupling(I, freqs[i])
       
        #Yig response
        fresp = inYig3(freqs, inpFFT, I, f0, ad)

        #outpYig = resp * coupling(I, f0)
        outpFFT = fresp
        for i in range(len(freqs)):
            outpFFT[i] = outpFFT[i] * coupling(I, freqs[i])

        #back to the time domain
        outpFFT = np.fft.ifftshift(outpFFT)
        outpYig = np.fft.ifft(outpFFT)

        #total output: Yig plus direct coupling
        output = outpYig + outpAntD
        outputr = output
        outputf = output

        if (I >= 2.1499 and I <=2.1509) or (I >= 2.2099 and I <= 2.2109):
            plt.plot(times, inp)
            plt.plot(times, inpAnt)
            plt.plot(times, np.fft.ifft(fresp))
            plt.show()
            #plt.plot(freqs, inpFFT0)
            plt.plot(freqs, inpFFT)
            plt.plot(freqs, fresp)
            #plt.plot(times, output)
            plt.show()
            #return
        
        #filtering
        lp = 3.5e9
        hp = 7.2e9
        #plt.plot(freqs, np.fft.fftshift(np.fft.fft(output)))
        #plt.show()
        b, a = signal.butter(5, [2*lp/sf, 2*hp/sf], btype = 'bandpass')
        output = signal.lfilter(b, a, output)

        #plotr(output)
        
        #mixing
        #plt.plot(freqs, np.fft.rfft(output))
        output *= np.sin(LO*2*np.pi*times)
        #plt.plot(freqs, np.fft.rfft(output))
        #plt.show()

        #Only filter around 
        df2 = 1e6
        b, a = signal.butter(12, [2*(f0 - df2)/sf, 2*(f0 + df2)/sf], btype = 'bandpass')
        outputf = signal.lfilter(b, a, outputf)


        #more filtering
        lp = 1.0e9
        b, a = signal.butter(5, 2*lp/sf, btype = 'lowpass')
        output = signal.lfilter(b, a, output)

        #sample the data: Voltage only!
        samp = 2.5e9
        mult = sf/samp
        #output = np.real(output) # voltage only
        relSamp = int(samp/sf*len(output))
        output = np.array([output[mult*i] for i in range(relSamp)])
        ntimes = np.array([times[mult*i] for i in range(relSamp)])

        #digital downconversion
        sf2 = 1./abs(ntimes[1]-ntimes[0])
        output *= np.exp(-1j*IF*2*np.pi*ntimes)
        #TODO this isnt entirely correct
        
        #more filtering
        lp = 62.5e6
        b, a = signal.butter(5, 2*lp/sf2, btype = 'lowpass')
        output = signal.lfilter(b, a, output)

        ##more filtering
        #lp = 62.5e6
        #b, a = signal.butter(5, 2*lp/sf2, btype = 'lowpass')
        #output = signal.lfilter(b, a, output)
        
        data.append(list(output))
        datar.append(list(outputr))
        dataf.append(list(outputf))


    fTime = 800
    data = np.array(data)
    pdata = data[:,:fTime]

    cmap = plt.cm.gnuplot
    times = np.arange(0, 800e-9, 0.4e-9)

    #raw
    im0 = plt.imshow(np.transpose(abs(np.real(datar))), aspect='auto', origin='lower', interpolation='none',cmap=cmap)
    im0.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    ax0 = im0.get_axes()
    ax0.set_title('Raw Data')
    plt.show()
    
    #bandpass filtered
    im02 = plt.imshow(np.transpose(abs(np.real(dataf))), aspect='auto', origin='lower', cmap=cmap, interpolation='none')
    im02.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    ax02 = im02.get_axes()
    ax02.set_title('bandPass filtered')
    plt.show()
    
    #logplot
    #im01 = plt.imshow(np.transpose(np.log(abs(np.real(datar)))), aspect='auto', origin='lower', cmap=cmap, interpolation='none')
    #ax01 = im01.get_axes()
    #ax01.set_title('logPlot raw')
    #im01.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    #plt.show()
    
    #filtered etc
    data = np.array(pdata)
    frdat = np.real(pdata)
    fidat = np.imag(pdata)
    fsqdat = np.sqrt(np.abs(frdat*np.conj(frdat)) + np.abs(fidat*np.conj(fidat)))
    fsqmdat = np.sqrt(abs(np.abs(frdat*np.conj(frdat)) - np.abs(fidat*np.conj(fidat))))

    #detail plot
    fig, (ax1, ax2) = plt.subplots(2)
    im03 = ax1.imshow(np.transpose(fsqdat), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im04 = ax2.imshow(np.transpose(np.log(fsqdat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    ax1.set_title('sqrt(I^2 + Q^1), filtered')
    ax2.set_title('real(dat), filtered')
    im03.set_extent([currents[0], currents[-1], ntimes[0], times[fTime]])
    im04.set_extent([currents[0], currents[-1], ntimes[0], times[fTime]])

 
    fig, ((ax1, ax2),(ax3, ax4), (ax5, ax6),(ax7, ax8)) = plt.subplots(4,2)
    im1 = ax1.imshow(np.transpose(frdat), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im2 = ax2.imshow(np.transpose(fidat), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im3 = ax3.imshow(np.transpose(fsqdat), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im4 = ax4.imshow(np.transpose(fsqmdat), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im5 = ax5.imshow(np.transpose(abs(frdat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im6 = ax6.imshow(np.transpose(abs(fidat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im6 = ax6.imshow(np.transpose(abs(fidat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im7 = ax7.imshow(np.transpose(np.log(fsqdat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im8 = ax8.imshow(np.transpose(np.log(fsqmdat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im1.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im2.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im3.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im4.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im5.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im6.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im7.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im8.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    ax1.set_title('real')
    ax2.set_title('imag')
    ax3.set_title('i^2 + q^2')
    ax4.set_title('i^2 - q^2')
    ax5.set_title('abs(real)')
    ax6.set_title('abs(imag)')
    ax7.set_title('log ( i^2 + q^2 )')
    ax8.set_title('log ( i^2 - q^2 )')
    plt.suptitle('Filtered and stuff')
    plt.show()

    #data
    rdat = np.real(datar)
    idat = np.imag(datar)
    sqdat = np.sqrt(np.abs(rdat*np.conj(rdat)) + np.abs(idat*np.conj(idat)))
    sqmdat = np.sqrt(np.abs(rdat*np.conj(rdat)) - np.abs(idat*np.conj(idat)))

    fig, ((ax1, ax2),(ax3, ax4), (ax5, ax6)) = plt.subplots(3,2)
    im1 = ax1.imshow(np.transpose(rdat), cmap=cmap, aspect=500000, origin='lower', interpolation='none')
    im2 = ax2.imshow(np.transpose(idat), cmap=cmap, aspect=500000, origin='lower', interpolation='none')
    im3 = ax3.imshow(np.transpose(sqdat), cmap=cmap, aspect=500000, origin='lower', interpolation='none')
    im4 = ax4.imshow(np.transpose(sqmdat), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im5 = ax5.imshow(np.transpose(abs(rdat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im6 = ax6.imshow(np.transpose(abs(idat)), cmap=cmap, aspect=500000, origin='lower',interpolation='none')
    im1.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im2.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im3.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im4.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im5.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    im6.set_extent([currents[0], currents[-1], ntimes[0], ntimes[fTime]])
    ax1.set_title('real')
    ax2.set_title('imag')
    ax3.set_title('i^2 + q^2')
    ax4.set_title('i^2 - q^2')
    ax5.set_title('abs(real)')
    ax6.set_title('abs(imag)')
    plt.suptitle('Before DC etc')
    plt.show()
    #return pdata
     
def inYig3(freqs, inputF, I, f0, d=ad):
    '''
    Tries to find response of YIG to full frequency spectrum of the input signal.
    '''
    
    tTime = 1./(freqs[1]-freqs[0])

    #determine the YIG passband (in terms of our downconverted frequencies)
    wH = wHI(I)
    wMin = np.sqrt(wH*(wH+wM)) 
    wMax = abs(wH + 0.5*wM)
    output = np.zeros(len(inputF))
    newWave = []

    #find the time delay due to the magnon propagation
    for i in range(len(freqs)):
        f = freqs[i]
        w = 2*np.pi*f
        #A time delay by phase-adding in the FFT is a valid method for shifting the ENTIRE time series, so the shift should be the same for the entire frequency series as well.
        #check the frequency is in the range of allowed MSSW frequencies
        xF = 0*2*pi*1e9
        #if ((w >= wMin and (w < wMax)) or ((w>-wMax) and (w<= -wMin))):
        #if ((w >= wMin and (w < wMax + xF)) or ((w>-wMax-xF) and (w<= -wMin))):
        #if (((w >= wMin - xF) and (w < wMax + xF)) or ((w>-wMax-xF) and (w<= -wMin +xF))):
        if abs(w) < wMax:
            td = pulseDelay(I, f)#/tTime
            #td = pulseDelay(I, f0, d) #TODO The above line caused the strange behaviour
            #Dont doubleCount the delay by using w here again!
            #TODO Warning! This is cast to real numbers!
            inputF[i] = np.complex(inputF[i]) * np.exp(-1j * f * td)
            #inputF[i] = inputF[i] * np.exp(1j * vPhase(I, f) * td)
            newWave.append(inputF[i])
        else:
            newWave.append(0)

    newWave = np.array(newWave)
    #before transforming back, shift fft back
    #newWave = np.fft.ifftshift(newWave)
    #newWave = np.fft.ifft(newWave)
    #plt.plot(ti, amps)
    #plt.plot(ti, newWave)
    #plt.show()
    
    return newWave
   

def averages(av=5):
    data = main(sin=True,noise=True)
    times = np.arange(0, 796e-9, .4e-9)
    I0 = 2.2
    currents = np.arange(I0-.15, I0+.25, .001)
    for i in range(av-1):
        data += main(sin=True,noise=True)

    cmap = plt.cm.gnuplot
    im = plt.imshow(np.transpose(data), cmap=cmap, aspect=500000, origin='lower')
    im.set_extent([currents[0], currents[-1], times[0], times[-1]])
    plt.show()

    #couple back out of yig

    #mix with output antenna

    #in the meantime: direct coupling, use data and interpolation?








#=====================================================================================================================
#---------------------------------------------PLOTTING----------------------------------------------------------------
#=====================================================================================================================


times = np.arange(0, 100e-9, 1e-9)

def plot2dim():
    plt.plot(times, pulse(times))
    plt.xlabel('time (s)')
    plt.show()

def plotr2(xdata, ydata, title='None'):
    plt.plot(xdata, ydata)
    plt.title(title)
    plt.show()

def plotr(xdata, title='None'):
    plt.plot(xdata)
    plt.title(title)
    plt.show()




def main2(pulse='squ'):
    '''
    Try to do the simulation more fully - simulate the YIG interaction with much higher time resolution
    and then apply all the measurement-like filtering that we actually use.
    '''

    #couplingConstant = .1
    
    times = np.arange(0, 800e-9, .04e-9)
    sf = 1./(times[1]-times[0])
    print 'sampling frequency', sf/1e9, ' GHz'
    freqs = np.linspace(-sf/2,sf/2.,len(times))
    #generate input pulse
    if pulse == 'squ':
        inp = squarePulseSin(times)
    elif pulse == 'dc':
        inp = squarePulse(times)
    elif pulse == 'gauss':
        inp = gaussPulse(times)

    #mix with input antenna
    atype = 'rlc'  #antenna type, can be rlc or log
    inpAnt = inputAnt(times, inp, BW=300e6, type=atype) #BW is ignored for rlc filters
    #plt.plot(inpAnt)
    #plt.show()
    #plt.plot(freqs, np.fft.rfft(inpAnt))
    #plt.show()
    
    #For direct coupling: Mix with box response and output antenna
    boxSig = sampleBoxResponse(times, inpAnt)
    outpAntD = inputAnt(times, boxSig, type=atype)

    plotr2(times, outpAntD)
    return
    #plot fft
    outDF = np.fft.fft(outpAntD)
    outDF = np.fft.fftshift(outDF)
    #print len(outDF)
    tf = np.arange(1, .5/(times[1]-times[0]) + 1, 1./times[-1])
    #print len(tf)
    #plotr(outDF)

    #couple into yig
    I0 = 2.2
    currents = np.arange(I0-.15, I0+.2, .01)
    data = []
    datar = []
    dataf = []
    #currents = [2.2,2.31]
    for I in currents:
        #input coupling due to antenna shape
        #inpAnt1 = coupling(I,f0) * inpAnt
        inpFFT = np.fft.fft(inpAnt)
        inpFFT = np.fft.fftshift(inpFFT)
        #inpFFT0 = inpFFT.copy()
        for i in range(len(freqs)):
            inpFFT[i] = inpFFT[i] * coupling(I, freqs[i])
       
        #Yig response
        fresp = inYig2(freqs, inpFFT, I, f0, ad)

        #outpYig = resp * coupling(I, f0)
        outpFFT = fresp
        for i in range(len(freqs)):
            outpFFT[i] = outpFFT[i] * coupling(I, freqs[i])

        #back to the time domain
        outpFFT = np.fft.ifftshift(outpFFT)
        outpYig = np.fft.ifft(outpFFT)

        #total output: Yig plus direct coupling
        output = outpYig + outpAntD
        outputr = output
        outputf = output

        if (I >= 2.1499 and I <=2.1509) or (I >= 2.2099 and I <= 2.2109):
            plt.plot(times, inp)
            plt.plot(times, inpAnt)
            plt.plot(times, np.fft.ifft(fresp))
            plt.show()
            #plt.plot(freqs, inpFFT0)
            plt.plot(freqs, inpFFT)
            plt.plot(freqs, fresp)
            #plt.plot(times, output)
            plt.show()
            #return
        
        #filtering
        lp = 3.5e9
        hp = 8.5e9
        #plt.plot(freqs, np.fft.fftshift(np.fft.fft(output)))
        #plt.show()
        b, a = signal.butter(5, [2*lp/sf, 2*hp/sf], btype = 'bandpass')
        output = signal.lfilter(b, a, output)

        #plotr(output)
        
        #mixing
        #plt.plot(freqs, np.fft.rfft(output))
        output *= np.sin(LO*2*np.pi*times)
        #plt.plot(freqs, np.fft.rfft(output))
        #plt.show()

        #Only filter around 
        df2 = 20
        b, a = signal.butter(12, [f0/sf - df2, f0/sf + df2], btype = 'bandpass')
        outputf = signal.lfilter(b, a, outputf)


        #more filtering
        lp = 1.0e9
        b, a = signal.butter(5, 2*lp/sf, btype = 'lowpass')
        output = signal.lfilter(b, a, output)

        #sample the data
        relSamp = int(2.5e9/sf*len(output))
        output = np.array([output[8*i] for i in range(relSamp)])

        #digital downconversion
        times2 = np.arange(0, 800e-9, 0.4e-9)
        sf2 = 1./(times2[1]-times2[0])
        output *= np.sin(IF*2*np.pi*times2)
        #TODO this isnt entirely correct
        
        #more filtering
        lp = 62.5e6
        b, a = signal.butter(5, 2*lp/sf2, btype = 'lowpass')
        output = signal.lfilter(b, a, output)
        
        data.append(list(abs(output)))
        datar.append(list(abs(outputr)))
        dataf.append(list(abs(outputf)))


    cmap = plt.cm.gnuplot
    times = np.arange(0, 800e-9, 0.4e-9)

    #raw
    im0 = plt.imshow(np.transpose(datar), aspect='auto', origin='lower', cmap=cmap)
    im0.set_extent([currents[0], currents[-1], times[0], times[800]])
    plt.show()
    
    #bandpass filtered
    im02 = plt.imshow(np.transpose(dataf), aspect='auto', origin='lower', cmap=cmap)
    im02.set_extent([currents[0], currents[-1], times[0], times[800]])
    ax02 = im02.get_axes()
    ax02.set_title('bandPass filtered')
    plt.show()
    
    #logplot
    im01 = plt.imshow(np.transpose(np.log(datar)), aspect='auto', origin='lower', cmap=cmap)
    ax01 = im01.get_axes()
    ax01.set_title('logPlot raw')
    im01.set_extent([currents[0], currents[-1], times[0], times[800]])
    plt.show()
    
    #plot a subset
    data = np.array(data)
    pdata = data[:,0:1600]
    im = plt.imshow(np.transpose(pdata), cmap=cmap, aspect=500000, origin='lower')
    im.set_extent([currents[0], currents[-1], times[0], times[800]])
    plt.show()
    return pdata
    

def inYig2(freqs, inputF, I, f0, d=ad):
    '''
    Tries to find response of YIG to full frequency spectrum of the input signal.
    '''
    
    #determine the YIG passband (in terms of our downconverted frequencies)
    wH = wHI(I)
    wMin = np.sqrt(wH*(wH+wM)) 
    wMax = abs(wH + 0.5*wM)
    output = np.zeros(len(inputF))
    newWave = []

    #find the time delay due to the magnon propagation
    for i in range(len(freqs)):
        w = 2*np.pi*freqs[i]
        #A time delay by phase-adding in the FFT is a valid method for shifting the ENTIRE time series, so the shift should be the same for the entire frequency series as well.
        td = pulseDelay(I, f0, d)
        #check the frequency is in the range of allowed MSSW frequencies
        if ((w >= wMin and (w < wMax)) or ((w>-wMax) and (w<= -wMin))):
            #td = pulseDelay(I, f0, d) #TODO The above line caused the strange behaviour
            #Dont doubleCount the delay by using w here again!
            inputF[i] = inputF[i] * np.exp(-1j * freqs[i] * td)
            inputF[i] = inputF[i] * np.exp(1j * vPhase(I, freqs[i]) * td)
            newWave.append(inputF[i])
        else:
            newWave.append(0)

    newWave = np.array(newWave)
    #before transforming back, shift fft back
    newWave = np.fft.ifftshift(newWave)
    newWave = np.fft.ifft(newWave)
    #plt.plot(ti, amps)
    #plt.plot(ti, newWave)
    #plt.show()
    
    return newWave
   


