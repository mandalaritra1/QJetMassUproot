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
        const = parameters[2]
        return mean,sigma,const
    
    def storeParameters(self):
        for i in range(len(self.histList)):
            hist = self.histList[i][0]
            parameters = self.fitGauss(hist)
            self.gaussMeanList[i] = parameters[0]
            self.gaussWidthList[i] = parameters[1]
            self.gaussConstList[i] = parameters[2]
            self.gaussMeanErrList[i] = parameters[3]
            parameters = self.fitbw(hist)
            self.bwMeanList[i] = parameters[0]
            self.bwWidthList[i] = parameters[1]
            self.bwConstList[i] = parameters[2]
            
    def plotGaussparameters(self):
        plt.figure(figsize = (16,5))
        plt.subplot(131)
        plt.errorbar(10*np.arange(self.N),self.gaussMeanList, self.gaussMeanErrList)
        plt.xlabel("pt,GeV")
        plt.ylabel("Mean")
        plt.subplot(132)
        plt.plot(10*np.arange(self.N),self.gaussWidthList)
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
        parameters, covariance = curve_fit(self.gauss, self.frac_values, self.histList[n][0],bounds = ([0.5,0.05,0.2],[2,0.3,20])) #fitting
        parameters1, covariance1 = curve_fit(self.bw, self.frac_values, self.histList[n][0]) #fitting_bw

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

        plt.plot(self.frac_values,self.histList[n][0], 'r',label = 'data')
        plt.legend()
    