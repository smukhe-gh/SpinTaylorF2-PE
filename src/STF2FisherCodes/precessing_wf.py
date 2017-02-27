import numpy as np
from numpy import sqrt, sin, cos, pi, exp, arctan2, arccos, array
from lal import MTSUN_SI
import pycbc
from pycbc.types import TimeSeries, FrequencySeries, zeros
from pycbc.waveform import get_fd_waveform, get_td_waveform
#import Cuda_envi
#from Cuda_envi import ctx_gpu, ctx_cpu

def m1m2_to_mchirpeta(mass1, mass2):
	eta = mass1*mass2/(mass1+mass2)/(mass1+mass2)
	mc = (mass1+mass2) * pow(eta, 3./5.)
	return (mc, eta)

def mchirpeta_to_m1m2(mchirp, eta):
	mtot = mchirp * pow(eta, -3./5)
	fac = sqrt(1. - 4.*eta)
	return (mtot * (1. + fac) / 2., mtot * (1. - fac) / 2.)

def rotateY(lst, angle):
	sinR, cosR = sin(angle), cos(angle)
	return [lst[0]*cosR + lst[2]*sinR, lst[1], -lst[0]*sinR + lst[2]*cosR]

def rotateZ(lst, angle):
	sinR, cosR = sin(angle), cos(angle)
	return [lst[0]*cosR - lst[1]*sinR, lst[0]*sinR + lst[1]*cosR, lst[2]]

def to_lal_coords(m1, m2, chi1, kappa, thetaJ, psiJ, alpha0, f0):
	""" Converts Andy and Richard's preferred coordinates to the antiquated LAL convention

	Returns inclination, psi0, S1hat vector"""
	v0 = pow(pi*MTSUN_SI*(m1+m2)*f0, 1./3.)
	gamma = m1*chi1*v0/m2
	denom = sqrt(1. + 2.*kappa*gamma + gamma*gamma)
	sinB, cosB = gamma*sqrt(1.-kappa*kappa)/denom, (1. + kappa*gamma)/denom
	sinA, cosA = sin(alpha0), cos(alpha0)

	lhat = [sinB*cosA, sinB*sinA, cosB]
	shat = [-sinB*cosA/gamma, -sinB*sinA/gamma, (kappa+gamma)/denom]

	lhat = rotateZ(rotateY(lhat, thetaJ), psiJ)
	shat = rotateZ(rotateY(shat, thetaJ), psiJ)

	psi0 = arctan2(lhat[1], lhat[0])
	incl = arccos(lhat[2])

	shat = rotateZ(shat, -psi0)
	shat = rotateY(shat, -incl) # Additional rotation present here.

	return (incl, psi0, array(shat))

class template:
	def __init__(self, m1, m2, chi1, kappa, thetaJ, psiJ, alpha0, phi0, f_inj):
		#m1, m2 = mchirpeta_to_m1m2(mchirp, eta)
		incl, psi0, shat = to_lal_coords(m1, m2, chi1, kappa, thetaJ, psiJ, alpha0, f_inj)

		self.mass1, self.mass2 = m1, m2
		self.spin1x, self.spin1y, self.spin1z = chi1*shat
		self.spin2x, self.spin2y, self.spin2z = (0., 0., 0.)
		self.distance = 1.
		self.inclination = incl
		self.pol = psi0
		self.coa_phase = phi0

class waveform:
	def __init__(self, f_inj, f_max, delta_f, amplitude_order=0, phase_order=7, spin_order=5, approximant='SpinTaylorF2', **kwargs):
		self._finj = f_inj
		self._deltaf = delta_f
		self._nsamples = int(f_max/delta_f)+1
		self._approximant = approximant
		self._ampO = amplitude_order
		self._phaseO = phase_order
		self._spinO = spin_order
		self._kwargs = kwargs
	def waveform(self, m1=2., m2=1., chi1=0., kappa=1., thetaJ=0.05, psiJ=0.05, alpha0=0., phi0=0., tC=0.):

		template_params = template(m1, m2, chi1, kappa, thetaJ, psiJ, alpha0, phi0, self._finj)
		hp, hx = get_fd_waveform(template_params,
		                         approximant=self._approximant,
		                         delta_f=self._deltaf,
		                         f_lower=self._finj,
		                         phase_order=self._phaseO,
		                         spin_order=self._spinO,
                                 sideband=None,
		                         amplitude_order=self._ampO,**(self._kwargs))
		#print('time elapsed=%f'%t.elapsed)

		sin2Y, cos2Y = sin(2.*template_params.pol), cos(2.*template_params.pol)
		wave = pycbc.DYN_RANGE_FAC* (hp*cos2Y+hx*sin2Y)
		wave.resize(self._nsamples)
		shift = FrequencySeries(exp(2.*pi*(0.+1.j)*tC*self._deltaf*np.arange(len(wave), dtype=np.complex128)), delta_f=self._deltaf, copy=False)

		return wave*shift

class fisher:
	def __init__(self, psd, wf_gen, f_low, f_high=None):
		self._psd = psd
		self._wfgen = wf_gen
		self._flow = f_low
		self._fhigh = f_high
	def _sigmasq(self, h1):
		return pycbc.filter.sigmasq(h1, psd=self._psd, low_frequency_cutoff=self._flow, high_frequency_cutoff=self._fhigh)
	def _overlap(self, h1, h2, normalized=False):
		return pycbc.filter.overlap(h1, h2, psd=self._psd, low_frequency_cutoff=self._flow, high_frequency_cutoff=self._fhigh, normalized=normalized)
	def _deriv(self, wf_params, key, dval):
		pcpy = wf_params.copy()
		pcpy[key] += dval
		h2 = self._wfgen.waveform(**pcpy)
		pcpy = wf_params.copy()
		pcpy[key] -= dval
		h1 = self._wfgen.waveform(**pcpy)
		return (h2-h1)/(2.*dval)
	def _calc_derivs(self, wf_params, wf_derivs, deriv_lst):
		h0 = self._wfgen.waveform(**wf_params)
		norm = 1./self._sigmasq(h0)
		derivs = {'norm': norm}
		for key in deriv_lst:
			dval = wf_derivs[key]
			derivs[key] = self._deriv(wf_params, key, dval)
		return derivs
	def calc_matrix(self, wf_params, wf_derivs, deriv_lst):
		derivs = self._calc_derivs(wf_params, wf_derivs, deriv_lst)
		result = np.zeros((len(deriv_lst), len(deriv_lst)))
		norm = derivs['norm']
		for ii, k1 in enumerate(deriv_lst):
			for jj, k2 in enumerate(deriv_lst):
				if ii <= jj:
					result[ii,jj] = norm*self._overlap(derivs[k1], derivs[k2])
				else:
					result[ii,jj] = result[jj,ii]
		return result
	def test_derivs(self, wf_params, wf_derivs, deriv_lst):
		derivs1 = self._calc_derivs(wf_params, wf_derivs, deriv_lst)
		half_wf_derivs = dict([(key, 0.5*val) for key, val in wf_derivs.items()])
		derivs2 = self._calc_derivs(wf_params, half_wf_derivs, deriv_lst)
		return dict([(key, self._overlap(derivs1[key], derivs2[key], normalized=True)) for key in deriv_lst])

