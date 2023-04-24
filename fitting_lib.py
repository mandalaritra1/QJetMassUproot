from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

class Histfit:
    def __init__(self, histList, frac_values, binlength):
        self.frac_values = frac_values
        self.binlength = binlength
        self.histList = histList
        self.N = len(histList)
        
        self.gaussMeanList = np.zeros(self.N)
        self.gaussWidthList = np.zeros(self.N)
        self.gaussConstList = np.zeros(self.N)
        self.gaussMeanErrList = np.zeros(self.N)
        
        self.bwMeanList = np.zeros(self.N)
        self.bwWidthList = np.zeros(self.N)
        self.bwConstList = np.zeros(self.N)
        
        self.gaussExpSigmaList = np.zeros(self.N)
        self.ptreco = [  15.,   25.,   35.,   45.,   55.,   65.,   75.,   85.,   95.,
        110.,  130.,  150.,  180., 200.]
        
        for i in range(len(self.histList)):
            self.histList[i] = histList[i]/(np.sum(histList[i])*self.binlength)
    def gauss(self,x,  x0, sigma,a):
        return (a*(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)))

    def bw(self,x, x0, gamma,k):
        """
        Returns the relativistic Breit-Wigner distribution for mass m,
        resonance mass m0 and width gamma.
        """
        #k = 2 * np.sqrt(2) * m0 * gamma * np.sqrt(m0**2 * (m0**2 + gamma**2))
        return k / ((x**2 - x0**2)**2 + (x**2)*gamma**2)
    
    def crystall_ball(self, x0, sigma, n, alpha, N):
        A = ((n / np.abs(alpha))**n ) * np.exp( - (np.abs(alpha)**2) / 2   )
        B = (n / np.abs(alpha)) - np.abs(alpha)
        #N = 1 / (sigma * (C+D))
        exp_part = N * np.exp (- ((x-x0)**2)/(2 * sigma**2) )
        polyn_part = N * A * (B - ((x - x0)/sigma) )**(-n)
        
        return np.where( ((x - x0)/sigma) > -alpha , exp_part, polyn_part)

    def gaussExp( self, x , x0 , sigma, k, N):
    
        A = N * np.exp( -0.5*( (x-x0)/sigma )**2)

        B = N * np.exp( ( k**2/2) - k*((x-x0)/sigma) )

        return np.where (  k < ((x-x0)/sigma)  , B, A)
    
    def fitGauss(self,hist):
        parameters, covariance = curve_fit(self.gauss, self.frac_values, hist,bounds = ([0.5,0.05,-5],[2,0.3,20]))
        mean = parameters[0]
        sigma = parameters[1]
        const = parameters[2]
        meanErr = covariance[0][0]
        return mean,sigma,const, meanErr
    
    def fitbw(self,hist):
        parameters, covariance = curve_fit(self.bw, self.frac_values, hist,bounds = ([0.5,0.05,-5],[2,0.3,5]))
        mean = parameters[0]
        sigma = parameters[1]
        k = parameters[2]
        return mean,sigma,k
    
    def fitGaussExp(self, hist):
        parameters, covariance = curve_fit(self.gaussExp, self.frac_values, hist)
        mean = parameters[0]
        sigma = parameters[1]
        const = parameters[3]
        k = parameters[2]
        return mean,sigma,k, const
    
    def storeParameters(self,skipFirst = False):
        if (skipFirst == True):
            for i in range(1,len(self.histList)):
                hist = self.histList[i]
                parameters = self.fitGauss(hist)
                self.gaussMeanList[i] = parameters[0]
                self.gaussWidthList[i] = parameters[1]
                self.gaussConstList[i] = parameters[2]
                self.gaussMeanErrList[i] = parameters[3]
                parameters = self.fitbw(hist)
                self.bwMeanList[i] = parameters[0]
                self.bwWidthList[i] = parameters[1]
                self.bwConstList[i] = parameters[2]

                parameters1 = self.fitGaussExp(hist)
                self.gaussExpSigmaList[i] = parameters1[1]
        else:
            for i in range(len(self.histList)):
                hist = self.histList[i]
                parameters = self.fitGauss(hist)
                self.gaussMeanList[i] = parameters[0]
                self.gaussWidthList[i] = parameters[1]
                self.gaussConstList[i] = parameters[2]
                self.gaussMeanErrList[i] = parameters[3]
                parameters = self.fitbw(hist)
                self.bwMeanList[i] = parameters[0]
                self.bwWidthList[i] = parameters[1]
                self.bwConstList[i] = parameters[2]

                parameters1 = self.fitGaussExp(hist)
                self.gaussExpSigmaList[i] = parameters1[1]

    def plotGaussExpParameters(self):
        plt.figure(figsize = (8,5))
        plt.plot(self.ptreco,self.gaussExpSigmaList,'bo-')
        plt.xlabel('pT (GeV)')
        plt.ylabel(r'$\sigma_{reco}/\sigma_{gen}$')
        
    def plotGaussparameters(self):
        plt.figure(figsize = (16,5))
        plt.subplot(131)
        plt.errorbar(10*np.arange(self.N),self.gaussMeanList, self.gaussMeanErrList)
        plt.xlabel("pt,GeV")
        plt.ylabel("Mean")
        plt.subplot(132)
        plt.plot(self.ptreco,self.gaussWidthList, 'bo-')
        plt.xlabel("pt,GeV")
        plt.ylabel("Width")
        plt.subplot(133)
        plt.plot(10*np.arange(self.N),self.gaussConstList)
        plt.xlabel("pt,GeV")
        plt.ylabel("Const")
    def plotBWparameters(self):
        plt.figure(figsize = (16,5))
        plt.subplot(131)
        plt.plot(10*np.arange(self.N),self.bwMeanList)
        plt.xlabel("pt,GeV")
        plt.ylabel("Mean")
        plt.subplot(132)
        plt.plot(10*np.arange(self.N),self.bwWidthList)
        plt.xlabel("pt,GeV")
        plt.ylabel("Width")
        plt.subplot(133)
        plt.plot(10*np.arange(self.N),self.bwConstList)
        plt.xlabel("pt,GeV")
        plt.ylabel("Const")
    

     
    def showFit(self,n):
        parameters, covariance = curve_fit(self.gauss, self.frac_values, self.histList[n],bounds = ([0.5,0.05,0.2],[2,0.3,20])) #fitting
        parameters1, covariance1 = curve_fit(self.bw, self.frac_values, self.histList[n]) #fitting_bw

        mean = parameters[0]
        sigma = parameters[1]
        const = parameters[2]
        print(const)

        plt.figure(figsize = (8,5))

        #fit_dist2 = self.bw(self.frac_values,parameters1[0],parameters1[1],parameters1[2])
        #plt.plot(self.frac_values,fit_dist2,'g--',label = 'bw fit')

        #plt.plot(self.frac_values,self.histList[n][0], 'r',label = 'data')
        #plt.legend()



        fit_dist = self.gauss(self.frac_values,mean,sigma,const)
        plt.plot(self.frac_values,fit_dist,'b--',label = 'gauss fit')

        plt.plot(self.frac_values,self.histList[n], 'r',label = 'data')
        plt.legend()
        
    def showFitTail(self, n):
        parameters, covariance = curve_fit(self.gaussExp, self.frac_values, self.histList[n]) #bounds = ([0.5,0.05,0.2,0.2],[2,0.3,20,20])
        mean = parameters[0]
        sigma = parameters[1]
        const = parameters[3]
        k = parameters[2]
        
        plt.figure(figsize = (8,5))    
        fit_dist = self.gaussExp(self.frac_values,mean,sigma,k, const)

        plt.plot(self.frac_values,fit_dist,'b--',label = 'gaussExp fit')
        
        #plt.plot(self.frac_values,fit_dist_gauss,'g--',label = 'gaussian core')

        plt.plot(self.frac_values,self.histList[n], 'r',label = 'data')
        plt.legend()
        
    def showFitTogether(self,ni = 0, nf = 12):
        fig = plt.figure(figsize = (20,10))
        ptreco = [10,20,30,40,50,60,70,80,90,100,120,140,160,200,7000]
        # loop over the rows and columns
        for i in range(3):
            for j in range(4):
                parameters, covariance = curve_fit(self.gaussExp, self.frac_values, self.histList[i*4+j+1]) #bounds = ([0.5,0.05,0.2,0.2],[2,0.3,20,20])
                mean = parameters[0]
                sigma = parameters[1]
                const = parameters[3]
                k = parameters[2]
  
                fit_dist = self.gaussExp(self.frac_values,mean,sigma,k, const)
                fit_dist_gauss = self.gauss(self.frac_values,mean,sigma,0.4)
                
                ax = plt.subplot(3, 4, i*4+j+1)
                n = i*4+j+1
                ax.text(0.0, 0.2, str(ptreco[n-1])+' < pT < ' +str(ptreco[n]), fontsize=12, va = 'bottom')

                ax.plot(self.frac_values,fit_dist,'b--',label = 'fit')
                ax.set_xlabel(r"$p_{T,jet}/p_{T,Z}$")
                ax.set_ylabel("Counts ( Normalized )")

                #plt.plot(self.frac_values,fit_dist_gauss,'g--',label = 'gaussian core')

                ax.plot(self.frac_values,self.histList[i*4+j+1], 'r',label = 'data')
                ax.legend()
                plt.savefig('Plots/fit_together.png',dpi = 300)
        
                # create a subplot at position i*6+j+1
                # plot some data
        # adjust the spacing between subplots
        plt.tight_layout()
        # show the figure
        plt.show()
