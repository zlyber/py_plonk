import gmpy2
import copy
import itertools
from domain import Radix2EvaluationDomain
from transcript import transcript
from composer import StandardComposer
from transcript import transcript
from plonk_core.lookup import multiset
from plonk_core.src.permutation import mod 
from plonk_core.src.proof_system.prover_key import Prover_Key
from plonk_core.src.proof_system.pi import into_dense_poly
from plonk_core.src.proof_system import quotient_poly
from plonk_core.src.proof_system import linearisation_poly
from field import field
from arithmetic import NTT,INTT,from_coeff_vec,resize
from load import read_scalar_data
from KZG import kzg10
from bls12_381 import fq,fr


def gen_proof(pp, pk: Prover_Key, cs: StandardComposer, transcript: transcript.Transcript):
    #init Fr params (FFTfield)
    Fr=fr.FrParameters()
    #get FFT domain
    
    domain=Radix2EvaluationDomain.new(cs.circuit_bound(),Fr)
    n=domain.size
    transcript.append_pi(b"pi")

    #1. Compute witness Polynomials
    w_l_scalar=read_scalar_data("w_l_scalar.txt")
    w_r_scalar=read_scalar_data("w_r_scalar.txt")
    w_o_scalar=read_scalar_data("w_o_scalar.txt")
    w_4_scalar=read_scalar_data("w_4_scalar.txt")

    # w_l_poly = INTT(domain,w_l_scalar)
    # w_r_poly = INTT(domain,w_r_scalar)
    # w_o_poly = INTT(domain,w_o_scalar)
    # w_4_poly = INTT(domain,w_4_scalar)

    w_l_poly = from_coeff_vec(INTT(domain,w_l_scalar))
    w_r_poly = from_coeff_vec(INTT(domain,w_r_scalar))
    w_o_poly = from_coeff_vec(INTT(domain,w_o_scalar))
    w_4_poly = from_coeff_vec(INTT(domain,w_4_scalar))

    w_polys = [kzg10.LabeledPoly.new(label="w_l_poly",hiding_bound=None,poly=w_l_poly),kzg10.LabeledPoly.new(label="w_r_poly",hiding_bound=None,poly=w_r_poly),
               kzg10.LabeledPoly.new(label="w_o_poly",hiding_bound=None,poly=w_o_poly),kzg10.LabeledPoly.new(label="w_4_poly",hiding_bound=None,poly=w_4_poly)]
    w_commits, w_rands = kzg10.commit_poly(pp,w_polys,Fr)

    transcript.append(b"w_l",w_commits[0].commitment.value)
    transcript.append(b"w_r",w_commits[1].commitment.value)
    transcript.append(b"w_o",w_commits[2].commitment.value)
    transcript.append(b"w_4",w_commits[3].commitment.value)
    #2. Derive lookup polynomials

    # Generate table compression factor
    zeta = transcript.challenge_scalar(b"zeta",Fr)
    transcript.append(b"zeta",zeta)

    # Compress lookup table into vector of single elements
    t_multiset = multiset.MultiSet([pk.lookup.table_1,pk.lookup.table_2,
                           pk.lookup.table_3,pk.lookup.table_4])
    compressed_t_multiset = t_multiset.compress(zeta)
    #Compute table poly
    compressed_t_poly = INTT(domain,compressed_t_multiset.elements)
    table_poly = from_coeff_vec(compressed_t_poly)

    # Compute query table f
    # When q_lookup[i] is zero the wire value is replaced with a dummy
    # value currently set as the first row of the public table
    # If q_lookup[i] is one the wire values are preserved
    # This ensures the ith element of the compressed query table
    # is an element of the compressed lookup table even when
    # q_lookup[i] is 0 so the lookup check will pass

    q_lookup_pad = [gmpy2.mpz(0)] * (n - len(cs.q_lookup))
    padded_q_lookup = cs.q_lookup + q_lookup_pad

    f_scalars = multiset.MultiSet([[],[],[],[]])
    for q_lookup, w_l, w_r, w_o, w_4 in zip(padded_q_lookup, w_l_scalar, w_r_scalar, w_o_scalar, w_4_scalar):
        if q_lookup == 0:
            f_scalars.elements[0].append(compressed_t_multiset.elements[0])
            for key in range(1,4):
                    f_scalars.elements[key].append(field(gmpy2.mpz(0),Fr))  
        else:
            f_scalars.elements[0].append(w_l)
            f_scalars.elements[1].append(w_r)
            f_scalars.elements[2].append(w_o)
            f_scalars.elements[3].append(w_4)

    # Compress all wires into a single vector
    compressed_f_multiset = f_scalars.compress(zeta)

    # Compute query poly
    compressed_f_poly = INTT(domain,compressed_f_multiset.elements)
    f_poly = from_coeff_vec(compressed_f_poly)
    f_polys = [kzg10.LabeledPoly.new(label="f_poly",hiding_bound=None,poly=f_poly)]

    # Commit to query polynomial
    f_poly_commit, _ = kzg10.commit_poly(pp,f_polys,Fr)
    transcript.append(b"f",f_poly_commit[0].commitment.value)

    # Compute s, as the sorted and concatenated version of f and t
    h_1, h_2 = compressed_t_multiset.combine_split(compressed_f_multiset)

    # Compute h polys
    h_1_temp = INTT(domain,h_1)
    h_2_temp = INTT(domain,h_2)
    h_1_poly = from_coeff_vec(h_1_temp)
    h_2_poly = from_coeff_vec(h_2_temp)

    # Commit to h polys
    h_1_polys = [kzg10.LabeledPoly.new(label="h_1_poly",hiding_bound=None,poly=h_1_poly)]
    h_2_polys = [kzg10.LabeledPoly.new(label="h_1_poly",hiding_bound=None,poly=h_2_poly)]
    h_1_poly_commit,_ = kzg10.commit_poly(pp,h_1_polys,Fr)
    h_2_poly_commit,_ = kzg10.commit_poly(pp,h_2_polys,Fr)

    # Add h polynomials to transcript
    transcript.append(b"h1", h_1_poly_commit[0].commitment.value)
    transcript.append(b"h2", h_2_poly_commit[0].commitment.value)

    # 3. Compute permutation polynomial

    # Compute permutation challenge `beta`.
    beta = transcript.challenge_scalar(b"beta",Fr)
    transcript.append(b"beta", beta)
    # Compute permutation challenge `gamma`.
    gamma = transcript.challenge_scalar(b"gamma",Fr)
    transcript.append(b"gamma", gamma)
    # Compute permutation challenge `delta`.
    delta = transcript.challenge_scalar(b"delta",Fr)
    transcript.append(b"delta", delta)
    # Compute permutation challenge `epsilon`.
    epsilon = transcript.challenge_scalar(b"epsilon",Fr)
    transcript.append(b"epsilon", epsilon)

    # Challenges must be different
    assert beta.value != gamma.value, "challenges must be different"
    assert beta.value != delta.value, "challenges must be different"
    assert beta.value != epsilon.value, "challenges must be different"
    assert gamma.value != delta.value, "challenges must be different"
    assert gamma.value != epsilon.value, "challenges must be different"
    assert delta.value != epsilon.value, "challenges must be different"
    
    z_poly = mod.compute_permutation_poly(domain,
        (w_l_scalar, w_r_scalar, w_o_scalar, w_4_scalar),
        beta,
        gamma,
        (
            pk.permutation.left_sigma[0],
            pk.permutation.right_sigma[0],
            pk.permutation.out_sigma[0],
            pk.permutation.fourth_sigma[0]
        ))
    # Commit to permutation polynomial.
    z_polys = [kzg10.LabeledPoly.new(label="z_poly",hiding_bound=None,poly=z_poly)]
    z_poly_commit,_ = kzg10.commit_poly(pp,z_polys,Fr)

    # Add permutation polynomial commitment to transcript.
    transcript.append(b"z", z_poly_commit[0].commitment.value)
    
    # Compute mega permutation polynomial.
    # Compute lookup permutation poly
    z_2_poly = mod.compute_lookup_permutation_poly(
        domain,
        compressed_f_multiset.elements,
        compressed_t_multiset.elements,
        h_1,
        h_2,
        delta,
        epsilon
    )

    # Commit to lookup permutation polynomial.
    z_2_polys = [kzg10.LabeledPoly.new(label="z_2_poly",hiding_bound=None,poly=z_2_poly)]
    z_2_poly_commit,_ = kzg10.commit_poly(pp,z_2_polys,Fr)

    # 3. Compute public inputs polynomial
    pi_poly = into_dense_poly(cs.public_inputs,cs.intended_pi_pos,n,Fr)

    # 4. Compute quotient polynomial

    # Compute quotient challenge `alpha`, and gate-specific separation challenges.
    alpha = transcript.challenge_scalar(b"alpha",Fr)
    transcript.append(b"alpha", alpha)

    range_sep_challenge = transcript.challenge_scalar(b"range separation challenge",Fr)
    transcript.append(b"range seperation challenge", range_sep_challenge)

    logic_sep_challenge = transcript.challenge_scalar(b"logic separation challenge",Fr)
    transcript.append(b"logic seperation challenge", logic_sep_challenge)

    fixed_base_sep_challenge = transcript.challenge_scalar(b"fixed base separation challenge",Fr)
    transcript.append(b"fixed base separation challenge", fixed_base_sep_challenge)

    var_base_sep_challenge = transcript.challenge_scalar(b"variable base separation challenge",Fr)
    transcript.append(b"variable base separation challenge", var_base_sep_challenge)

    lookup_sep_challenge = transcript.challenge_scalar(b"lookup separation challenge",Fr)
    transcript.append(b"lookup separation challenge", lookup_sep_challenge)

    t_poly = quotient_poly.compute(
        domain,pk,
        z_poly,z_2_poly,
        w_l_poly,w_r_poly,w_o_poly,w_4_poly,
        pi_poly,
        f_poly,table_poly,h_1_poly,h_2_poly,
        alpha,beta,gamma,delta,epsilon,zeta,
        range_sep_challenge,logic_sep_challenge,
        fixed_base_sep_challenge,
        var_base_sep_challenge,
        lookup_sep_challenge)

    t_i_poly = split_tx_poly(n, t_poly, Fr)

    t_i_polys = [kzg10.LabeledPoly.new(label="t_i_polys[0]",hiding_bound=None,poly=t_i_poly[0]),
               kzg10.LabeledPoly.new(label="t_i_polys[1]",hiding_bound=None,poly=t_i_poly[1]),
               kzg10.LabeledPoly.new(label="t_i_polys[2]",hiding_bound=None,poly=t_i_poly[2]),
               kzg10.LabeledPoly.new(label="t_i_polys[3]",hiding_bound=None,poly=t_i_poly[3]),
               kzg10.LabeledPoly.new(label="t_i_polys[4]",hiding_bound=None,poly=t_i_poly[4]),
               kzg10.LabeledPoly.new(label="t_i_polys[5]",hiding_bound=None,poly=t_i_poly[5]),
               kzg10.LabeledPoly.new(label="t_i_polys[6]",hiding_bound=None,poly=t_i_poly[6]),
               kzg10.LabeledPoly.new(label="t_i_polys[7]",hiding_bound=None,poly=t_i_poly[7])]
    
    t_commits, _ = kzg10.commit_poly(pp,t_i_polys,Fr)

    # Add quotient polynomial commitments to transcript
    transcript.append(b"t_1", t_commits[0].commitment.value)
    transcript.append(b"t_2", t_commits[1].commitment.value)
    transcript.append(b"t_3", t_commits[2].commitment.value)
    transcript.append(b"t_4", t_commits[3].commitment.value)
    transcript.append(b"t_5", t_commits[4].commitment.value)
    transcript.append(b"t_6", t_commits[5].commitment.value)
    transcript.append(b"t_7", t_commits[6].commitment.value)
    transcript.append(b"t_8", t_commits[7].commitment.value)

    # 4. Compute linearisation polynomial

    # Compute evaluation challenge `z`.
    z_challenge = transcript.challenge_scalar(b"z", Fr)
    transcript.append(b"z", z_challenge)

    lin_poly, evaluations = linearisation_poly.compute(
            domain,
            pk,
            alpha,beta,gamma,delta,epsilon,zeta,
            range_sep_challenge,
            logic_sep_challenge,
            fixed_base_sep_challenge,
            var_base_sep_challenge,
            lookup_sep_challenge,
            z_challenge,
            w_l_poly,w_r_poly,w_o_poly,w_4_poly,
            t_i_poly[0],
            t_i_poly[1],
            t_i_poly[2],
            t_i_poly[3],
            t_i_poly[4],
            t_i_poly[5],
            t_i_poly[6],
            t_i_poly[7],
            z_poly,
            z_2_poly,
            f_poly,
            h_1_poly,
            h_2_poly,
            table_poly)
    
    # Add evaluations to transcript.
    # First wire evals
    transcript.append(b"a_eval", evaluations.wire_evals.a_eval)
    transcript.append(b"b_eval", evaluations.wire_evals.b_eval)
    transcript.append(b"c_eval", evaluations.wire_evals.c_eval)
    transcript.append(b"d_eval", evaluations.wire_evals.d_eval)

    # Second permutation evals
    transcript.append(b"left_sig_eval", evaluations.perm_evals.left_sigma_eval)
    transcript.append(b"right_sig_eval",evaluations.perm_evals.right_sigma_eval)
    transcript.append(b"out_sig_eval", evaluations.perm_evals.out_sigma_eval)
    transcript.append(b"perm_eval", evaluations.perm_evals.permutation_eval)

    # Third lookup evals
    transcript.append(b"f_eval", evaluations.lookup_evals.f_eval)
    transcript.append(b"q_lookup_eval", evaluations.lookup_evals.q_lookup_eval)
    transcript.append(b"lookup_perm_eval",evaluations.lookup_evals.z2_next_eval)
    transcript.append(b"h_1_eval", evaluations.lookup_evals.h1_eval)
    transcript.append(b"h_1_next_eval", evaluations.lookup_evals.h1_next_eval)
    transcript.append(b"h_2_eval", evaluations.lookup_evals.h2_eval)

    # Fourth, all evals needed for custom gates
    for label, eval in evaluations.custom_evals.vals:
        static_label = label.encode('utf-8')
        transcript.append(static_label, eval)

    # 5. Compute Openings using KZG10
    #
    # We merge the quotient polynomial using the `z_challenge` so the SRS
    # is linear in the circuit size `n`

    # Compute aggregate witness to polynomials evaluated at the evaluation
    # challenge `z`
    aw_challenge = transcript.challenge_scalar(b"aggregate_witness", Fr)

    # XXX: The quotient polynomials is used here and then in the
    # opening poly. It is being left in for now but it may not
    # be necessary. Warrants further investigation.
    # Ditto with the out_sigma poly.
    aw_polys = [kzg10.LabeledPoly.new(label="lin_poly",hiding_bound=None,poly=lin_poly),
               kzg10.LabeledPoly.new(label="prover_key.permutation.left_sigma.0.clone()",hiding_bound=None,poly=pk.permutation.left_sigma[0]),
               kzg10.LabeledPoly.new(label="prover_key.permutation.right_sigma.0.clone()",hiding_bound=None,poly=pk.permutation.right_sigma[0]),
               kzg10.LabeledPoly.new(label="prover_key.permutation.out_sigma.0.clone()",hiding_bound=None,poly=pk.permutation.out_sigma[0]),
               kzg10.LabeledPoly.new(label="f_poly",hiding_bound=None,poly=f_poly),
               kzg10.LabeledPoly.new(label="h_2_poly",hiding_bound=None,poly=h_2_poly),
               kzg10.LabeledPoly.new(label="table_poly",hiding_bound=None,poly=table_poly)]
    
    aw_commits, aw_rands = kzg10.commit_poly(pp,aw_polys,Fr)
    aw_opening = kzg10.open(
        pp,
        itertools.chain(aw_polys, w_polys),
        itertools.chain(aw_commits, w_commits),
        z_challenge,
        aw_challenge,
        itertools.chain(aw_rands, w_rands),
        None
    )

    saw_challenge = transcript.challenge_scalar(b"aggregate_witness", Fr)
    saw_polys = [kzg10.LabeledPoly.new(label="z_poly",hiding_bound=None,poly=z_poly),
               kzg10.LabeledPoly.new(label="w_l_poly",hiding_bound=None,poly=w_l_poly),
               kzg10.LabeledPoly.new(label="w_r_poly",hiding_bound=None,poly=w_r_poly),
               kzg10.LabeledPoly.new(label="w_4_poly",hiding_bound=None,poly=w_4_poly),
               kzg10.LabeledPoly.new(label="h_1_poly",hiding_bound=None,poly=h_1_poly),
               kzg10.LabeledPoly.new(label="z_2_poly",hiding_bound=None,poly=z_2_poly),
               kzg10.LabeledPoly.new(label="table_poly",hiding_bound=None,poly=table_poly)]
    
    saw_commits, saw_rands = kzg10.commit_poly(pp,saw_polys,Fr)
    saw_opening = kzg10.open(
        pp,
        saw_polys,
        saw_commits,
        z_challenge.mul(domain.element(1)),
        saw_challenge,
        saw_rands,
        None
    )

    Proof = kzg10.Proof(
            a_comm = w_commits[0].commitment.value,
            b_comm = w_commits[1].commitment.value,
            c_comm = w_commits[2].commitment.value,
            d_comm = w_commits[3].commitment.value,
            z_comm = saw_commits[0].commitment.value,
            f_comm = f_poly_commit[0].commitment.value,
            h_1_comm = h_1_poly_commit[0].commitment.value,
            h_2_comm = h_2_poly_commit[0].commitment.value,
            z_2_comm = z_2_poly_commit[0].commitment.value,
            t_1_comm = t_commits[0].commitment.value,
            t_2_comm = t_commits[1].commitment.value,
            t_3_comm = t_commits[2].commitment.value,
            t_4_comm = t_commits[3].commitment.value,
            t_5_comm = t_commits[4].commitment.value,
            t_6_comm = t_commits[5].commitment.value,
            t_7_comm = t_commits[6].commitment.value,
            t_8_comm = t_commits[7].commitment.value,
            aw_opening = aw_opening,
            saw_opening = saw_opening,
            evaluations = evaluations)
    return Proof


def split_tx_poly(n, t_x, params):
    buf:list = t_x[:]
    buf = resize(buf, n << 3, field.zero(params))
    return [
        from_coeff_vec(buf[0:n]),
        from_coeff_vec(buf[n:2 * n]),
        from_coeff_vec(buf[2 * n:3 * n]),
        from_coeff_vec(buf[3 * n:4 * n]),
        from_coeff_vec(buf[4 * n:5 * n]),
        from_coeff_vec(buf[5 * n:6 * n]),
        from_coeff_vec(buf[6 * n:7 * n]),
        from_coeff_vec(buf[7 * n:])
    ]


