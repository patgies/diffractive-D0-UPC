import numpy as np

# One-loop running of alpha_s from alpha_s(mZ), with Nf = 5 above the b-quark threshold mb and
# Nf = 4 below it (alpha_s continuous at mb):
# 1/alpha_s(mu) = 1/alpha_s(mu0) + 2*b0(Nf)*ln(mu/mu0), b0(Nf) = (33 - 2*Nf)/(12*pi)

alphas_mZ = 0.118
mZ        = 91.2
mb        = 4.8
Nc        = 3

def b0(nf):
    return (33 - 2*nf) / (12*np.pi)

def alphas_run(mu):
    inv_mb = 1.0/alphas_mZ + 2*b0(5)*np.log(mb/mZ)          # 1/alpha_s(mb), Nf = 5 from mZ
    inv = np.where(mu >= mb,
                   1.0/alphas_mZ + 2*b0(5)*np.log(mu/mZ),
                   inv_mb + 2*b0(4)*np.log(mu/mb))
    return 1.0 / inv

if __name__ == '__main__':
    print('One-loop running alpha_s  (Nf=5 above mb=%.1f GeV, 4 below; alpha_s(mZ)=%.3f)' % (mb, alphas_mZ))
    print('=' * 36)
    print('   mu [GeV]    alpha_s(mu)')
    print('-' * 36)
    for mu in [1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 91.2]:
        print('  %8.1f      %.4f' % (mu, alphas_run(mu)))
    print()

    mc = 1.5   
    pt_values = np.arange(0.1, 10.1, 0.2)

    print('alpha_s at mu = sqrt(pT^2 + mc^2),  mc=%.1f GeV' % mc)
    print('=' * 40)
    print('   pT [GeV]   mu [GeV]   alpha_s(mu)')
    print('-' * 40)
    for pt in pt_values:
        mu = np.sqrt(pt**2 + mc**2)
        print('  %8.1f   %8.4f    %.4f' % (pt, mu, alphas_run(mu)))
