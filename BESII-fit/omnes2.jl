using Interpolations, CSV, DataFrames; using LinearAlgebra, StaticArrays

# Coupled-channel Omnes matrix form
omnes_BernOrsay_df = DataFrame(CSV.File("OmnesMatrix_Johanna.csv"))
omnes_DaiPennington_df = DataFrame(CSV.File("OmnesMatrix_Stefan.dat"))
# tmatrix_DaiPennington_df = DataFrame(CSV.File("../omnes/TMatrix_Stefan.csv"))
# delta20_df = DataFrame(CSV.File("../omnes/delta20_10GeV.csv"))
# absomnes20_df = DataFrame(CSV.File("../omnes/AbsOmnes20_10GeV.csv"))
# omnes20_df = DataFrame(CSV.File("../omnes/Omnes20_eetoJpsipipi_mpipi_4260.csv"))

const ReOmnes11_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ReOmega_11))
const ImOmnes11_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ImOmega_11))
const ReOmnes12_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ReOmega_12))
const ImOmnes12_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ImOmega_12))
const ReOmnes21_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ReOmega_21))
const ImOmnes21_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ImOmega_21))
const ReOmnes22_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ReOmega_22))
const ImOmnes22_BO = LinearInterpolation(collect(omnes_BernOrsay_df.sqrt_s), collect(omnes_BernOrsay_df.ImOmega_22))

const ReOmnes11_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ReOmega_11))
const ImOmnes11_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ImOmega_11))
const ReOmnes12_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ReOmega_12))
const ImOmnes12_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ImOmega_12))
const ReOmnes21_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ReOmega_21))
const ImOmnes21_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ImOmega_21))
const ReOmnes22_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ReOmega_22))
const ImOmnes22_DP = LinearInterpolation(collect(omnes_DaiPennington_df.sqrt_s), collect(omnes_DaiPennington_df.ImOmega_22))
	
# const ReT11_DP = LinearInterpolation(collect(tmatrix_DaiPennington_df.sqrt_s), collect(tmatrix_DaiPennington_df.ReT11))
# const ImT11_DP = LinearInterpolation(collect(tmatrix_DaiPennington_df.sqrt_s), collect(tmatrix_DaiPennington_df.ImT11))
# const ReT12_DP = LinearInterpolation(collect(tmatrix_DaiPennington_df.sqrt_s), collect(tmatrix_DaiPennington_df.ReT12))
# const ImT12_DP = LinearInterpolation(collect(tmatrix_DaiPennington_df.sqrt_s), collect(tmatrix_DaiPennington_df.ImT12))
# const ReT22_DP = LinearInterpolation(collect(tmatrix_DaiPennington_df.sqrt_s), collect(tmatrix_DaiPennington_df.ReT22))
# const ImT22_DP = LinearInterpolation(collect(tmatrix_DaiPennington_df.sqrt_s), collect(tmatrix_DaiPennington_df.ImT22))

const cmatrix2 = zeros(ComplexF64, 2, 2);


"""
	omnes_cc(x; which = :DP)

Coupled-channel Omnes function. Channel-1: ππ; channel-2: KKbar.
* `x` is the ππ invariant mass in units of GeV.
* `which = :BO` or `:DP` corresponds to the Bern-Orsay and Dai-Pennington parameterizations, respectively; default to `:DP`
"""
function omnes_cc(x; which = :BO)
	if which == :BO
		cmatrix2[1,1] = ReOmnes11_BO(x) + 1im* ImOmnes11_BO(x)
		cmatrix2[2,1] = ReOmnes21_BO(x) + 1im* ImOmnes21_BO(x)
		cmatrix2[1,2] = ReOmnes12_BO(x) + 1im* ImOmnes12_BO(x)
		cmatrix2[2,2] = ReOmnes22_BO(x) + 1im* ImOmnes22_BO(x)
	elseif which == :DP
		cmatrix2[1,1] = ReOmnes11_DP(x) + 1im* ImOmnes11_DP(x)
		cmatrix2[2,1] = ReOmnes21_DP(x) + 1im* ImOmnes21_DP(x)
		cmatrix2[1,2] = ReOmnes12_DP(x) + 1im* ImOmnes12_DP(x)
		cmatrix2[2,2] = ReOmnes22_DP(x) + 1im* ImOmnes22_DP(x)
    else
        error_cc("`omnes(x; which = $which)` is not defined.")
	end
	return SMatrix{2,2}(cmatrix2)
end


"""
	tmatrix_cc(x)

Coupled-channel T-matrix using the Dai-Pennington parameterization. Channel-1: ππ; channel-2: KKbar.
* `x` is the c.m. energy in units of GeV.
"""
# function tmatrix_cc(x)
# 	cmatrix2[1,1] = ReT11_DP(x) + 1im* ImT11_DP(x)
# 	cmatrix2[2,1] = ReT12_DP(x) + 1im* ImT12_DP(x)
# 	cmatrix2[1,2] = cmatrix2[2,1]
# 	cmatrix2[2,2] = ReT22_DP(x) + 1im* ImT22_DP(x)
# 	return SMatrix{2,2}(cmatrix2)
# end

#1204.6251
const omnes11re = DataFrame!(CSV.File("Omnes-1204-6251/Re11-1204.csv"));
const omnes11im = DataFrame!(CSV.File("Omnes-1204-6251/Im11-1204.csv"));
const omnes12re = DataFrame!(CSV.File("Omnes-1204-6251/Re12-1204.csv"));
const omnes12im = DataFrame!(CSV.File("Omnes-1204-6251/Im12-1204.csv"));
const omnes21re = DataFrame!(CSV.File("Omnes-1204-6251/Re21-1204.csv"));
const omnes21im = DataFrame!(CSV.File("Omnes-1204-6251/Im21-1204.csv"));
const omnes22re = DataFrame!(CSV.File("Omnes-1204-6251/Re22-1204.csv"));
const omnes22im = DataFrame!(CSV.File("Omnes-1204-6251/Im22-1204.csv"));

const omnes11intre=LinearInterpolation(omnes11re.x,omnes11re.y)
const omnes11intim=LinearInterpolation(omnes11im.x,omnes11im.y)
const omnes12intre=LinearInterpolation(omnes12re.x,omnes12re.y)
const omnes12intim=LinearInterpolation(omnes12im.x,omnes12im.y)
const omnes21intre=LinearInterpolation(omnes21re.x,omnes21re.y)
const omnes21intim=LinearInterpolation(omnes21im.x,omnes21im.y)
const omnes22intre=LinearInterpolation(omnes22re.x,omnes22re.y)
const omnes22intim=LinearInterpolation(omnes22im.x,omnes22im.y)

omnes11fun(w)=omnes11intre(w)+im*omnes11intim(w)
omnes12fun(w)=omnes12intre(w)+im*omnes12intim(w)
omnes21fun(w)=omnes21intre(w)+im*omnes21intim(w)
omnes22fun(w)=omnes22intre(w)+im*omnes22intim(w)


# Master thesis Johanna
const omnesjohannare = DataFrame!(CSV.File("Omnes-Johanna/Re.csv"));
const omnesjohannaim = DataFrame!(CSV.File("Omnes-Johanna/Im.csv"));

const omnesjohannaintre=LinearInterpolation(omnesjohannare.x,omnesjohannare.y)
const omnesjohannaintim=LinearInterpolation(omnesjohannaim.x,omnesjohannaim.y)

omnesjohannadata(w)=omnesjohannaintre(w)+im*omnesjohannaintim(w)

nothing
