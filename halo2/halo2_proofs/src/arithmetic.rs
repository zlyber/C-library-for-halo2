//! This module provides common utilities, traits and structures for group,
//! field and polynomial arithmetic.

//use std::intrinsics::size_of;
use std::ptr;
use super::multicore;
pub use ff::Field;

use lazy_static::lazy_static;
use std::fs::File;
use std::io::Write;
use std::sync::Mutex;
use std::time::Instant; //timer
use group::{
    ff::{BatchInvert, PrimeField},
    Curve, Group as _,
};

pub use halo2curves::{CurveAffine, CurveExt, FieldExt, Group};

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));


fn multiexp_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C], acc: &mut C::Curve) {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();
    let c = if bases.len() < 4 {
        1
    } else if bases.len() < 32 {
        3
    } else {
        (f64::from(bases.len() as u32)).ln().ceil() as usize
    };

    fn get_at<F: PrimeField>(segment: usize, c: usize, bytes: &F::Repr) -> usize {
        let skip_bits = segment * c;     
        let skip_bytes = skip_bits / 8;  
        
        
        if skip_bytes >= 32 {
            return 0;
        }

        let mut v = [0; 8];
        for (v, o) in v.iter_mut().zip(bytes.as_ref()[skip_bytes..].iter()) {
            *v = *o;
        }

        let mut tmp = u64::from_le_bytes(v);
        tmp >>= skip_bits - (skip_bytes * 8);
        tmp = tmp % (1 << c);

        tmp as usize
    }

    let segments = (256 / c) + 1;

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            *acc = acc.double();
        }

        #[derive(Clone, Copy)]
        enum Bucket<C: CurveAffine> {
            None,
            Affine(C),
            Projective(C::Curve),
        }

        impl<C: CurveAffine> Bucket<C> {
            fn add_assign(&mut self, other: &C) {
                *self = match *self {
                    Bucket::None => Bucket::Affine(*other),
                    Bucket::Affine(a) => Bucket::Projective(a + *other),
                    Bucket::Projective(mut a) => {
                        a += *other;
                        Bucket::Projective(a)
                    }
                }
            }

            fn add(self, mut other: C::Curve) -> C::Curve {
                match self {
                    Bucket::None => other,
                    Bucket::Affine(a) => {
                        other += a;
                        other
                    }
                    Bucket::Projective(a) => other + &a,
                }
            }
        }

        let mut buckets: Vec<Bucket<C>> = vec![Bucket::None; (1 << c) - 1];

        for (coeff, base) in coeffs.iter().zip(bases.iter()) {
            let coeff = get_at::<C::Scalar>(current_segment, c, coeff);
            if coeff != 0 {
                buckets[coeff - 1].add_assign(base);
            }
        }
        
        // Summation by parts
        // e.g. 3a + 2b + 1c = a +
        //                    (a) + b +
        //                    ((a) + b) + c

        let mut running_sum = C::Curve::identity();
        for exp in buckets.into_iter().rev() {
            running_sum = exp.add(running_sum);
            *acc = *acc + &running_sum;
        }
    }
}

/// Performs a small multi-exponentiation operation.
/// Uses the double-and-add algorithm with doublings shared across points.
pub fn small_multiexp<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();
    let mut acc = C::Curve::identity();

    // for byte idx
    for byte_idx in (0..32).rev() {
        // for bit idx
        for bit_idx in (0..8).rev() {
            acc = acc.double();
            // for each coeff
            for coeff_idx in 0..coeffs.len() {
                let byte = coeffs[coeff_idx].as_ref()[byte_idx];
                if ((byte >> bit_idx) & 1) != 0 {
                    acc += bases[coeff_idx];
                }
            }
        }
    }

    acc
}

/// Performs a multi-exponentiation operation.(MSM)
///
/// This function will panic if coeffs and bases have a different length.
///
/// This will use multithreading if beneficial.
pub fn best_multiexp<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
    
    let mut basex: Vec<Vec<u64>> = vec![vec![0; 4]; bases.len()];
    let mut basey: Vec<Vec<u64>> = vec![vec![0; 4]; bases.len()];
    let mut basez: Vec<Vec<u64>> = vec![vec![0; 4]; bases.len()];
    let mut Coeffs: Vec<Vec<u64>> = vec![vec![0; 4]; bases.len()];
    let coeff_len:u64=bases.len().try_into().unwrap();
    let mut acc_x:Vec<u64>=vec![0;4];
    let mut acc_y:Vec<u64>=vec![0;4];
    let mut acc_z:Vec<u64>=vec![0;4];
    let mut bases_mut = bases.to_vec();//额外的内存分配和数据复制，但目前没有想到避免的办法
    let mut coeffs_mut=coeffs.to_vec();
    unsafe{
    let mut acc = C::Curve::identity();
    let (prefix, aligned, suffix) = bases_mut.align_to_mut::<u64>();
    let (prefix, caligned, suffix) = coeffs_mut.align_to_mut::<u64>();
    for (i, chunk) in caligned.chunks_exact(4).enumerate() {
        Coeffs[i] = chunk.to_vec();
    }
    for (i, chunk) in aligned.chunks_exact(8).enumerate() {
        let mut chunk_iter = chunk.chunks_exact(4);
        basex[i] = chunk_iter.next().unwrap().to_vec();
        basey[i] = chunk_iter.next().unwrap().to_vec();
    }
    let mut coeff_ptr: Vec<*mut u64> = Coeffs.iter_mut()
              .map(|row| row.as_mut_ptr())
              .collect();
    let mut basex_ptr: Vec<*mut u64> = basex.iter_mut()
            .map(|row| row.as_mut_ptr())
            .collect();
    let mut basey_ptr: Vec<*mut u64> = basey.iter_mut()
        .map(|row| row.as_mut_ptr())
        .collect();
    let mut basez_ptr: Vec<*mut u64> = basez.iter_mut()
            .map(|row| row.as_mut_ptr())
            .collect();
    let accx_ptr = acc_x.as_mut_ptr();
    let accy_ptr = acc_y.as_mut_ptr();
    let accz_ptr = acc_z.as_mut_ptr();
    MSM(coeff_ptr.as_mut_ptr(), basex_ptr.as_mut_ptr(), basey_ptr.as_mut_ptr(), basez_ptr.as_mut_ptr(), accx_ptr,accy_ptr,accz_ptr,coeff_len);
    let mut acc_c: Vec<u64> = Vec::new();
    acc_c.extend_from_slice(&acc_x);
    acc_c.extend_from_slice(&acc_y);
    acc_c.extend_from_slice(&acc_z);
    let (prefix, ACC, suffix) = acc_c.align_to_mut::<C::Curve>();
    acc=ACC[0];
    
    acc
}
}

/// Performs a radix-$2$ Fast-Fourier Transformation (FFT) on a vector of size
/// $n = 2^k$, when provided `log_n` = $k$ and an element of multiplicative
/// order $n$ called `omega` ($\omega$). The result is that the vector `a`, when
/// interpreted as the coefficients of a polynomial of degree $n - 1$, is
/// transformed into the evaluations of this polynomial at each of the $n$
/// distinct powers of $\omega$. This transformation is invertible by providing
/// $\omega^{-1}$ in place of $\omega$ and dividing each resulting field element
/// by $n$.
///
/// This will use multithreading if beneficial.
pub fn best_fft<G: Group>(a: &mut [G], omega: G::Scalar, log_n: u32) {
    
    let p=omega.get_lower_128();

    let mut pomega: [u64; 4] = [0; 4];
    //K=5
    //omega_inv
    if p==318428984803038777672377917081039968157u128{
        pomega[0]=2898883370368880564;
        pomega[1]=426384372937323306;
        pomega[2]=14682641521443947722;
        pomega[3]=3462788798779817102;
    }
    //extended_omega
    else if p==82380007769639746197868146856976870763u128 {
        pomega[0]=12327492402416538496;
        pomega[1]=10101357140220402583;
        pomega[2]=8672023191981693263;
        pomega[3]=2067983392542172724;
    }
    //extended_omega_inv
    else {
        pomega[0]=2179530681542097973;
        pomega[1]=691179081718012817;
        pomega[2]=13745101439070792779;
        pomega[3]=2314654573598063593;
    }
    // //omega
    //     pomega[0]=2455860039615204452;
    //     pomega[1]=2952429719635721714;
    //     pomega[2]=16534267467014174157;
    //     pomega[3]=3135063178989759747;
    // 
    //K=10
    // //omega_inv
    // if p==297649304259851281220676526985611293155u128{
    //     pomega[0]=559035126084790881;
    //     pomega[1]=12123427087958374242;
    //     pomega[2]=2949720136027289278;
    //     pomega[3]=298260759310561979;
    // }
    // //extended_omega
    // else if p==75322075251223052345800817392097324774u128 {
    //     pomega[0]=11791621636447142361;
    //     pomega[1]=3213422488462342693;
    //     pomega[2]=7137044954843475233;
    //     pomega[3]=1120461048903492910;
    // }
    // //extended_omega_inv
    // else {
    //     pomega[0]=18046099818638051508;
    //     pomega[1]=17490710434920251023;
    //     pomega[2]=5069583672902197219;
    //     pomega[3]=2226268784200321758;
    // }
    // //omega
    //     pomega[0]=11338367819501839635;
    //     pomega[1]=3740034585862428569;
    //     pomega[2]=9914360941292543340;
    //     pomega[3]=2496945838366253767;
    //K=15
    // //omega_inv
    // if p==81462561015711620142297117044210472484u128{
    //     pomega[0]=16940376884292359863;
    //     pomega[1]=2219943494317580003;
    //     pomega[2]=948315165668358246;
    //     pomega[3]=2847967631620404566;
    // }
    // //extended_omega
    // else if p==308268835236386324909053348514485158890u128 {
    //     pomega[0]=10328996839774307296;
    //     pomega[1]=3137371740708723579;
    //     pomega[2]=14560387799605567379;
    //     pomega[3]=2967229345595982336;
    // }
    // //extended_omega_inv
    // else if p==95044254358722968092221589431928008955u128 {
    //     pomega[0]=10576018011517079560;
    //     pomega[1]=9716248366216130094;
    //     pomega[2]=7306589758931143380;
    //     pomega[3]=3472459677849482257;
    // }
    // //omega
    //     pomega[0]=7883942902451040980;
    //     pomega[1]=3486698467654550028;
    //     pomega[2]=6285631179351301021;
    //     pomega[3]=2186570735001675543;

    //K=20
    // //omega_inv
    // if p==237453209598177823829231244369545407839u128{
    //     pomega[0]=9091523228542722555;
    //     pomega[1]=13041570835191812554;
    //     pomega[2]=1589937418470923934;
    //     pomega[3]=1254658277821569380;
    // }
    // //extended_omega
    // else if p==340026751871900385633170862117082283742u128 {
    //     pomega[0]=5038868593009920868;
    //     pomega[1]=9032737560981761680;
    //     pomega[2]=13782952129176230699;
    //     pomega[3]=1728816712498300990;
    // }
    // //extended_omega_inv
    // else {
    //     pomega[0]=3881617269061619234;
    //     pomega[1]=2485871750576600637;
    //     pomega[2]=12653098317601961217;
    //     pomega[3]=604174078516457111;
    // }
    // //omega
    //     pomega[0]=138401318583551795;
    //     pomega[1]=9905347880094660640;
    //     pomega[2]=11718360747020793053;
    //     pomega[3]=1852013080597357056;
    // 
    //K=21
    // //omega_inv
    // if p==314774160241541424962438210521828737784u128{
    //     pomega[0]=7419436655583739328;
    //     pomega[1]=4546061169868866527;
    //     pomega[2]=17077186631929237990;
    //     pomega[3]=2356479378029159210;
    // }
    // //extended_omega
    // else if p==218895413970211164727092521504155764835u128 {
    //     pomega[0]=12738377379681167599;
    //     pomega[1]=2197391538012506479;
    //     pomega[2]=5373537405419956753;
    //     pomega[3]=1461168961025718367;
    // }
    // //extended_omega_inv
    // else {
    //     pomega[0]=3566522453088237343;
    //     pomega[1]=11152815455894748557;
    //     pomega[2]=11291125319747055124;
    //     pomega[3]=807261542184412802;
    // }
    // //omega
    //     pomega[0]=7063725442162364797;
    //     pomega[1]=10456422038513750871;
    //     pomega[2]=610002596162631866;
    //     pomega[3]=1686867808588488897;
    //K=22
    // //omega_inv
    // if p==299874562726571307606600617828231559524u128{
    //     pomega[0]=3881617269061619234;
    //     pomega[1]=2485871750576600637;
    //     pomega[2]=12653098317601961217;
    //     pomega[3]=604174078516457111;
    // }
    // //extended_omega
    // else if p==4565888064856207896893974357169208011u128 {
    //     pomega[0]=17104276958738319052;
    //     pomega[1]=7189382688254064745;
    //     pomega[2]=1604580299050760569;
    //     pomega[3]=2533986721453659612;
    // }
    // //extended_omega_inv
    // else {
    //     pomega[0]=6821884587516740863;
    //     pomega[1]=7826069202380984337;
    //     pomega[2]=14111528798780703700;
    //     pomega[3]=1014596263072848297;
    // }
    // //omega
    //     pomega[0]=5038868593009920868;
    //     pomega[1]=9032737560981761680;
    //     pomega[2]=13782952129176230699;
    //     pomega[3]=1728816712498300990;
    //K=23
    // //omega_inv
    // if p==262401740751341575011060224449330837291u128{
    //     pomega[0]=3566522453088237343;
    //     pomega[1]=11152815455894748557;
    //     pomega[2]=11291125319747055124;
    //     pomega[3]=807261542184412802;
    // }
    // //extended_omega
    // else if p==50737075762603355302678067946645737824u128 {
    //     pomega[0]=13338605924273364442;
    //     pomega[1]=11440449704248451096;
    //     pomega[2]=16859609365912477452;
    //     pomega[3]=3421252324365184758;
    // }
    // //extended_omega_inv
    // else {
    //     pomega[0]=2738242980467392064;
    //     pomega[1]=8765460162850139420;
    //     pomega[2]=6637814084492473216;
    //     pomega[3]=1260493707339115276;
    // }
    // //omega
    //     pomega[0]=12738377379681167599;
    //     pomega[1]=2197391538012506479;
    //     pomega[2]=5373537405419956753;
    //     pomega[3]=1461168961025718367;
    //K=24
    // //omega_inv
    // if p==320929653853306149548133271619207731414u128{
    //     pomega[0]=6821884587516740863;
    //     pomega[1]=7826069202380984337;
    //     pomega[2]=14111528798780703700;
    //     pomega[3]=1014596263072848297;
    // }
    // //extended_omega
    // else if p==48150040116460294936019791923596749975u128 {
    //     pomega[0]=8521028017761076590;
    //     pomega[1]=16651763990358389415;
    //     pomega[2]=8564895374847835975;
    //     pomega[3]=2458108212177418459;
    // }
    // //extended_omega_inv
    // else {
    //     pomega[0]=2612153849575289297;
    //     pomega[1]=10023378403827285182;
    //     pomega[2]=16006860937169638985;
    //     pomega[3]=118158886462801266;
    // }
    // //omega
    //     pomega[0]=17104276958738319052;
    //     pomega[1]=7189382688254064745;
    //     pomega[2]=1604580299050760569;
    //     pomega[3]=2533986721453659612;
    // 
    let omega_ptr = &mut pomega as *mut u64;
    unsafe{
    let (prefix, aligned, suffix) = a.align_to_mut::<[u64;4]>();//可行的方法

    NTT(aligned.as_mut_ptr() , omega_ptr , log_n);

    let (prefix, a, suffix) = aligned.align_to_mut::<G>();
    }

}

/// This perform recursive butterfly arithmetic
pub fn recursive_butterfly_arithmetic<G: Group>(
    a: &mut [G],
    n: usize,
    twiddle_chunk: usize,
    twiddles: &[G::Scalar],
) {
    if n == 2 {
        let t = a[1];
        a[1] = a[0];
        a[0].group_add(&t);
        a[1].group_sub(&t);
    } else {
        let (left, right) = a.split_at_mut(n / 2);
        rayon::join(
            || recursive_butterfly_arithmetic(left, n / 2, twiddle_chunk * 2, twiddles),
            || recursive_butterfly_arithmetic(right, n / 2, twiddle_chunk * 2, twiddles),
        );

        // case when twiddle factor is one
        let (a, left) = left.split_at_mut(1);
        let (b, right) = right.split_at_mut(1);
        let t = b[0];
        b[0] = a[0];
        a[0].group_add(&t);
        b[0].group_sub(&t);

        left.iter_mut()
            .zip(right.iter_mut())
            .enumerate()
            .for_each(|(i, (a, b))| {
                let mut t = *b;
                t.group_scale(&twiddles[(i + 1) * twiddle_chunk]);
                *b = *a;
                a.group_add(&t);
                b.group_sub(&t);
            });
    }
}

/// Convert coefficient bases group elements to lagrange basis by inverse FFT.
pub fn g_to_lagrange<C: CurveAffine>(g_projective: Vec<C::Curve>, k: u32) -> Vec<C> {
    let n_inv = C::Scalar::TWO_INV.pow_vartime(&[k as u64, 0, 0, 0]);
    let mut omega_inv = C::Scalar::ROOT_OF_UNITY_INV;
    for _ in k..C::Scalar::S {
        omega_inv = omega_inv.square();
    }

    let mut g_lagrange_projective = g_projective;
    best_fft(&mut g_lagrange_projective, omega_inv, k);
    parallelize(&mut g_lagrange_projective, |g, _| {
        for g in g.iter_mut() {
            *g *= n_inv;
        }
    });

    let mut g_lagrange = vec![C::identity(); 1 << k];
    parallelize(&mut g_lagrange, |g_lagrange, starts| {
        C::Curve::batch_normalize(
            &g_lagrange_projective[starts..(starts + g_lagrange.len())],
            g_lagrange,
        );
    });

    g_lagrange
}

/// This evaluates a provided polynomial (in coefficient form) at `point`.
pub fn eval_polynomial<F: Field>(poly: &[F], point: F) -> F {
    fn evaluate<F: Field>(poly: &[F], point: F) -> F {
        poly.iter()
            .rev()
            .fold(F::zero(), |acc, coeff| acc * point + coeff)
    }
    let n = poly.len();
    let num_threads = multicore::current_num_threads();
    if n * 2 < num_threads {
        evaluate(poly, point)
    } else {
        let chunk_size = (n + num_threads - 1) / num_threads;
        let mut parts = vec![F::zero(); num_threads];
        multicore::scope(|scope| {
            for (chunk_idx, (out, poly)) in
                parts.chunks_mut(1).zip(poly.chunks(chunk_size)).enumerate()
            {
                scope.spawn(move |_| {
                    let start = chunk_idx * chunk_size;
                    out[0] = evaluate(poly, point) * point.pow_vartime(&[start as u64, 0, 0, 0]);
                });
            }
        });
        parts.iter().fold(F::zero(), |acc, coeff| acc + coeff)
    }
}

/// This computes the inner product of two vectors `a` and `b`.
///
/// This function will panic if the two vectors are not the same size.
pub fn compute_inner_product<F: Field>(a: &[F], b: &[F]) -> F {
    // TODO: parallelize?
    assert_eq!(a.len(), b.len());

    let mut acc = F::zero();
    for (a, b) in a.iter().zip(b.iter()) {
        acc += (*a) * (*b);
    }

    acc
}

/// Divides polynomial `a` in `X` by `X - b` with
/// no remainder.
pub fn kate_division<'a, F: Field, I: IntoIterator<Item = &'a F>>(a: I, mut b: F) -> Vec<F>
where
    I::IntoIter: DoubleEndedIterator + ExactSizeIterator,
{
    b = -b;
    let a = a.into_iter();

    let mut q = vec![F::zero(); a.len() - 1];

    let mut tmp = F::zero();
    for (q, r) in q.iter_mut().rev().zip(a.rev()) {
        let mut lead_coeff = *r;
        lead_coeff.sub_assign(&tmp);
        *q = lead_coeff;
        tmp = lead_coeff;
        tmp.mul_assign(&b);
    }

    q
}

/// This simple utility function will parallelize an operation that is to be
/// performed over a mutable slice.
pub fn parallelize<T: Send, F: Fn(&mut [T], usize) + Send + Sync + Clone>(v: &mut [T], f: F) {
    let n = v.len();
    let num_threads = multicore::current_num_threads();
    let mut chunk = (n as usize) / num_threads;
    if chunk < num_threads {
        chunk = 1;
    }

    multicore::scope(|scope| {
        for (chunk_num, v) in v.chunks_mut(chunk).enumerate() {
            let f = f.clone();
            scope.spawn(move |_| {
                let start = chunk_num * chunk;
                f(v, start);
            });
        }
    });
}

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow + 1)) <= num {
        pow += 1;
    }

    pow
}

/// Returns coefficients of an n - 1 degree polynomial given a set of n points
/// and their evaluations. This function will panic if two values in `points`
/// are the same.
pub fn lagrange_interpolate<F: FieldExt>(points: &[F], evals: &[F]) -> Vec<F> {
    assert_eq!(points.len(), evals.len());
    if points.len() == 1 {
        // Constant polynomial
        vec![evals[0]]
    } else {
        let mut denoms = Vec::with_capacity(points.len());
        for (j, x_j) in points.iter().enumerate() {
            let mut denom = Vec::with_capacity(points.len() - 1);
            for x_k in points
                .iter()
                .enumerate()
                .filter(|&(k, _)| k != j)
                .map(|a| a.1)
            {
                denom.push(*x_j - x_k);
            }
            denoms.push(denom);
        }
        // Compute (x_j - x_k)^(-1) for each j != i
        denoms.iter_mut().flat_map(|v| v.iter_mut()).batch_invert();

        let mut final_poly = vec![F::zero(); points.len()];
        for (j, (denoms, eval)) in denoms.into_iter().zip(evals.iter()).enumerate() {
            let mut tmp: Vec<F> = Vec::with_capacity(points.len());
            let mut product = Vec::with_capacity(points.len() - 1);
            tmp.push(F::one());
            for (x_k, denom) in points
                .iter()
                .enumerate()
                .filter(|&(k, _)| k != j)
                .map(|a| a.1)
                .zip(denoms.into_iter())
            {
                product.resize(tmp.len() + 1, F::zero());
                for ((a, b), product) in tmp
                    .iter()
                    .chain(std::iter::once(&F::zero()))
                    .zip(std::iter::once(&F::zero()).chain(tmp.iter()))
                    .zip(product.iter_mut())
                {
                    *product = *a * (-denom * x_k) + *b * denom;
                }
                std::mem::swap(&mut tmp, &mut product);
            }
            assert_eq!(tmp.len(), points.len());
            assert_eq!(product.len(), points.len() - 1);
            for (final_coeff, interpolation_coeff) in final_poly.iter_mut().zip(tmp.into_iter()) {
                *final_coeff += interpolation_coeff * eval;
            }
        }
        final_poly
    }
}

pub(crate) fn evaluate_vanishing_polynomial<F: FieldExt>(roots: &[F], z: F) -> F {
    fn evaluate<F: FieldExt>(roots: &[F], z: F) -> F {
        roots.iter().fold(F::one(), |acc, point| (z - point) * acc)
    }
    let n = roots.len();
    let num_threads = multicore::current_num_threads();
    if n * 2 < num_threads {
        evaluate(roots, z)
    } else {
        let chunk_size = (n + num_threads - 1) / num_threads;
        let mut parts = vec![F::one(); num_threads];
        multicore::scope(|scope| {
            for (out, roots) in parts.chunks_mut(1).zip(roots.chunks(chunk_size)) {
                scope.spawn(move |_| out[0] = evaluate(roots, z));
            }
        });
        parts.iter().fold(F::one(), |acc, part| acc * part)
    }
}

pub(crate) fn powers<F: FieldExt>(base: F) -> impl Iterator<Item = F> {
    std::iter::successors(Some(F::one()), move |power| Some(base * power))
}

#[cfg(test)]
use rand_core::OsRng;

#[cfg(test)]
use crate::halo2curves::pasta::Fp;

#[test]
fn test_lagrange_interpolate() {
    let rng = OsRng;

    let points = (0..5).map(|_| Fp::random(rng)).collect::<Vec<_>>();
    let evals = (0..5).map(|_| Fp::random(rng)).collect::<Vec<_>>();

    for coeffs in 0..5 {
        let points = &points[0..coeffs];
        let evals = &evals[0..coeffs];

        let poly = lagrange_interpolate(points, evals);
        assert_eq!(poly.len(), points.len());

        for (point, eval) in points.iter().zip(evals) {
            assert_eq!(eval_polynomial(&poly, *point), *eval);
        }
    }
}
