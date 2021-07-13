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

nothing
