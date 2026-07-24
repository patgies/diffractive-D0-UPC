import numpy as np

# One-loop running of alpha_s
# Input: alpha_s(mZ) from PDG world average
# Formula: 1/alpha_s(mu) = 1/alpha_s(mZ) + 2*b0*ln(mu/mZ)
# b0 = (33 - 2*Nf) / (12*pi)  for Nc=3

alphas_mZ = 0.118    # PDG world average
mZ        = 91.2     # GeV
Nc        = 3
Nf        = 4        # active flavours at charm scale (u, d, s, c)
b0        = (33 - 2*Nf) / (12*np.pi)

def alphas_run(mu):
    return 1.0 / (1.0/alphas_mZ + 2*b0*np.log(mu/mZ))

if __name__ == '__main__':
    print('One-loop running alpha_s  (Nf=%d, alpha_s(mZ)=%.3f)' % (Nf, alphas_mZ))
    print('=' * 36)
    print('   mu [GeV]    alpha_s(mu)')
    print('-' * 36)
    for mu in [1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 91.2]:
        print('  %8.1f      %.4f' % (mu, alphas_run(mu)))
    print()

    # Print at transverse mass scale mu = sqrt(pT^2 + mc^2) 
    mc = 1.5   # charm mass GeV
    pt_values = np.arange(0.1, 10.1, 0.2)

    print('alpha_s at mu = sqrt(pT^2 + mc^2),  mc=%.1f GeV' % mc)
    print('=' * 40)
    print('   pT [GeV]   mu [GeV]   alpha_s(mu)')
    print('-' * 40)
    for pt in pt_values:
        mu = np.sqrt(pt**2 + mc**2)
        print('  %8.1f   %8.4f    %.4f' % (pt, mu, alphas_run(mu)))
