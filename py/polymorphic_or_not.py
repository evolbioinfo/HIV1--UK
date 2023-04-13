import pandas as pd

NONPOLYMORPHIC = 'nonpolymorphic'
POLYMORPHIC = 'polymorphic'
NNRTI_STANFORD_TEXT = """
V90I
V90I is a polymorphic accessory NNRTI-resistance mutation which occurs in 0.5% to 5.0% of ART-naïve persons depending on subtype. It is selected in vitro by EFV (1), ETR (1), and RPV (2). It is commonly selected in persons receiving RPV (3,4,5) and weakly selected in persons receiving NVP and EFV (4). Alone it causes little, if any, reduction in susceptibility to any NNRTI (6,7).

A98G
A98G is a relatively nonpolymorphic accessory NNRTI-resistance mutation. It occurs in 0.1% to 0.5% of untreated persons depending on subtype. It has been selected primarily in persons receiving NVP and EFV. It reduces NVP, EFV, RPV and DOR susceptibility about 2-fold (6,8,9,10,11).

L100I
L100I is a nonpolymorphic mutation selected in vitro by EFV (12,13), ETR (14), and RPV (2). It is selected in persons receiving EFV (15), ETR (16), and RPV (3). L100I rarely occurs in isolation. In combination with K103N, it causes >50-fold reduced susceptibility to NVP and EFV (6), 10-fold reduced susceptibility to RPV, and 5 to 10-fold reduced susceptibility to ETR (2,6,14,17) and DOR (11,18).

L100V is an extremely rare nonpolymorphic mutation selected in vitro by RPV (13). Its impact on NNRTI susceptibility has been assessed in few isolates (7).

K101E/H/P/Q/R/N
K101E is a nonpolymorphic mutation commonly selected in persons receiving each of the NNRTIs except possibly DOR for which few data are available (3,15,19,20). K101E usually occurs in combination with other NNRTI-resistance mutations. Alone it confers 3-10-fold reduced susceptibility to NVP and about 2-fold reduced EFV, ETR, and RPV susceptibility (2,3,6,13,18).

K101P is a nonpolymorphic two base-pair mutation selected in persons receiving each of the NNRTIs except DOR for which few data are available (3,15,19,20). Alone it reduces susceptibility to NVP, EFV, and RPV >20-fold and ETR susceptibility about 5-fold (2,6,19,21,22). It does not reduce DOR susceptibility (23,24).

K101H is a nonpolymorphic mutation weakly selected in persons receiving NVP, EFV, and ETR (19,25). When it occurs in combination with other NNRTI-resistance mutations, it contributes reduced susceptibility to NVP and EFV (6).

K101Q is a polymorphic accessory mutation occurring in 1% to 3% ART-naïve persons depending on subtype, It is weakly selected in patients receiving NVP and EFV (25). It appears to have minimal, if any, detectable effect on NNRTI susceptibility (6). K101N/A/T are nonpolymorphic NNRTI-selected mutations (25). Their effects on NNRTI susceptibility have not been studied. K101R is a polymorphic mutation that is not selected by NNRTIs. It has no effect on NNRTI susceptibility.

K103N/S/H/T/R/Q/E
K103N is a nonpolymorphic NNRTI-resistance mutation selected in persons receiving NVP and EFV (15,26,27,28). It reduces NVP and EFV susceptibility about 50 and 20-fold, respectively (6,29,30) but does not reduce susceptibility to RPV, ETR, or DOR (6,20). K103N is the most commonly transmitted drug-resistance mutation (31,32). The presence of K103N prior to starting the WHO first-line regimen tenofovir/lamivudine/efavirenz is associated with an increased risk of VF (33,34), although it appears to be only minimally elevated in those patients who have K103N and no other resistance mutation (35).

K103S is a nonpolymorphic two base-pair mutation selected by NVP and EFV. It usually occurs in samples from persons who previously had K103N (36). It is associated with about 5-fold reduced susceptibility to EFV and >20-fold reduced susceptibility to NVP (6,37) Because K103S is a 2-bp change from the wildtype K, persons with K103S may be more likely to harbor K103N, which is just a 1-bp change from wildtype. Indeed, a high proportion of viruses with K103S contain it as part of a mixture with K103N.

K103H is a rare nonpolymorphic NNRTI-selected mutation that confers about 20-fold reduced susceptibility to NVP and EFV (37,38). K103T is a rare nonpolymorphic NNRTI-selected mutation that confers >20-fold reduced susceptibility to NVP but does not appear to reduce susceptibility to EFV (9,37,38).

K103R is a polymorphic mutation that occurs in about 1% of ART-naïve persons. Alone, it has no effect on NNRTI susceptibility. However, in combination with V179D (and possibly V179E), K103R reduces NVP and EFV susceptibility about 15-fold (21). K103E/Q are rare mutations that do not appear to be selected by NNRTIs or to be associated with reduced susceptibility to any NNRTI.

V106A/M/I
V106A is a nonpolymorphic mutation selected in persons receiving NVP (36,39,40,41) and DOR (13,20,42). It confers about 50-fold reduced susceptibility to NVP and about 5-fold reduced susceptibility to EFV (6,9,38,43). Alone it causes about 10-fold reduced DOR susceptibility (10,18). In combination with other DOR-associated DRMs including mutations at positions 225, 227, and 234, it is associated with about 100-fold reduced DOR susceptibility (11,20,44).

V106M is a nonpolymorphic mutation selected primarily by EFV and NVP. It occurs more often in subtype C viruses because it requires a single base-pair change in subtype C viruses – GTG (V) => ATG (M) – but a two base-pair change in all other subtypes – GTA (V) => ATG (M) (45,46,47,48,49,50,51). It confers >30-fold reduced susceptibility to NVP and EFV (6,9). It is selected in vitro and in vivo by DOR and data on a single isolate suggests that alone it is associated with about 3-fold reduced DOR susceptibility (13,18,42).

V106I is a polymorphic accessory NNRTI-selected mutation that occurs in 0.1% to 3% of viruses from NNRTI-naïve persons depending on subtype. It has been reported to be selected commonly in persons receiving NVP and EFV and in 4 of 10 persons with VF and HIVDR while receiving DOR (20,42,52). It has been reported infrequently in persons receiving ETR and RPV (4,19). Alone it appears to have little if any effect on NNRTI susceptibility (10,53). In combination with V179D, it causes low-level reductions in susceptibility to NVP, EFV, ETR, and probably RPV (54). In combination with mutations at position 227, it is associated with high-levels of reduced susceptibility to DOR (20). In a weighted genotypic score, it was associated with a reduced virological response to treatment with an ETR-containing regimen when combined with other NNRTI-resistance mutations (8).

V108I
V108I is a relatively nonpolymorphic accessory mutation that is selected in vitro by NVP, EFV, and DOR (12,13,55). It has commonly been reported in persons receiving NVP, EFV, and ETR (15,19,36) and in one of 10 persons with VF and HIVDR in persons receiving DOR (42). It reduces NVP and EFV susceptibility about 2-fold but does not appear to reduce ETR, RPV or DOR susceptibility (6,10,53,56).

E138K/A/G/Q/R
E138K is a nonpolymorphic mutation that is selected in a high proportion of persons receiving RPV (3) and in a smaller proportion of persons receiving ETR (19). Alone it reduces RPV and ETR susceptibility about 2-fold; in combination with M184I or K101E it reduces RPV and ETR susceptibility about 3 fold (57,58,59). The combination of E138K and M184I or K101E appears sufficient to cause VF on an RPV-containing regimen (3). E138K does not appear to reduce susceptibility to NVP, EFV, or DOR (6,10).

E138A is a polymorphic mutation that ranges in prevalence from about 2% to 5% in ART-naïve persons depending on subtype. It is not known whether it is selected in NNRTI-treated persons. It reduces ETR and RPV susceptibility about 2-fold (2,6,18,57,60). In a weighted genotypic score, it was associated with a reduced virological response to treatment with an ETR-containing regimen when combined with other NNRTI-resistance mutations (8). Based on limited data, baseline E138 mutations do not appear to be associated with an increased risk of VF in persons receiving an RPV-containing regimen (61,62).

E138Q/G are nonpolymorphic accessory mutations selected by ETR (3,19,57) and occasionally by NVP and EFV. E138Q/G are associated with about 3-fold reduced NVP susceptibility and 2-fold reduced ETR and RPV susceptibility (6,57). E138R is an extremely rare 2 base-pair mutation that is selected primarily in non-subtype B viruses from NNRTI-treated persons. It appears to reduce NNRTI susceptibility slightly more than E138Q/G (2).

V179D/E/F/I/L/T
V179D is a polymorphic accessory mutation occurring in about 1% of ART-naïve persons. It is occasionally selected in persons receiving EFV. It reduces NVP, EFV, ETR, and RPV susceptibility about 2-fold but does not appear to reduce DOR susceptibility (6,18,38). In a weighted genotypic score, it was associated with a reduced virological response to treatment with an ETR-containing regimen when combined with other NNRTI-resistance mutations (8). The combination of V179D and K103R act synergistically to reduce NVP and EFV susceptibility about 15-fold (21). The combination of V179D and V106I also appear to act synergistically to reduce NVP and EFV susceptibility (54). The presence of V179D does not appear to interfere with the response to a first-line EFV containing regimen (63,64).

V179F is a nonpolymorphic accessory mutation that usually emerges in combination with Y181C. It is one of the most common mutations emerging in persons receiving ETR (19,65). Alone it does not reduce NNRTI susceptibility, but the combination of V179F and Y181C confers >10-fold reduced susceptibility to RPV and ETR (2,8).

V179E is a nonpolymorphic mutation weakly selected by EFV. Its effect on NNRTI susceptibility is not well defined (9,38). V179T is an uncommon nonpolymorphic mutation that is infrequently selected in persons receiving NNRTIs. Alone it does not appear to be associated with reduced NNRTI susceptibility. In a weighted genotypic score, it was associated with a reduced virological response to treatment with an ETR-containing regimen when combined with other NNRTI-resistance mutations (8).

V179L is a rare NNRTI-selected mutation that is listed as an RPV-associated drug resistance mutation in the RPV package insert. Its effect on NNRTI susceptibility has not been studied.

V179I is a highly polymorphic mutation occurring in 3% to 50% of viruses from ART-naïve persons depending on subtype. It is selected by each of the NNRTIs except possibly DOR (3,19,65).

Y181C/I/V
Y181C is a nonpolymorphic mutation selected in vitro by each of the NNRTIs except DOR (2,12,14,55,66). It is selected in vivo primarily by NVP, ETR, and RPV (3,19,36). It confers >50-fold reduced susceptibility to NVP (67,68), about 5-fold reduced susceptibility to ETR (8,38,69), about 3-fold reduced susceptibility to RPV (3), and <2-fold reduced susceptibility to EFV (6,68). Although Y181C alone minimally reduces EFV susceptibility, it has been associated with a reduced response to treatment with an EFV-containing regimen in persons with VF on an NVP-containing regimen presumably because other NNRTI-resistance mutations are often present as minority NNRTI-resistance variants. It has been associated with a reduced virological response to treatment with an ETR-containing regimen when combined with other NNRTI-resistance mutations (8). Y181C alone does not reduce DOR susceptibility (18).

Y181I/V are 2-base pair nonpolymorphic mutations selected by NVP and ETR (19). Y181I/V confer >50-fold reduced susceptibility to NVP and 0 to 15-fold reduced susceptibility to ETR and RPV (2). Y181I/V have been associated with a reduced virological response to treatment with an ETR-containing regimen (8). In clinical isolates, Y181V has been associated with variable reductions in DOR susceptibility (18).

Y181F/S/G are rare nonpolymorphic NNRTI-associated mutations that are usually present as part of an electrophoretic mixture. They likely represent transitional mutations between Y and I or V or possibly partial revertant mutations. Y181S reduces NVP susceptibility about 30-fold (30).

Y188L/C/H
Y188L is a nonpolymorphic 2-base pair mutation selected by NVP and EFV (15,36). It confers >50-fold reduced susceptibility to NVP and EFV (29,56,68,70), 5-fold reduced susceptibility to RPV (6,71) and minimally reduced susceptibility to ETR (6,72,73). Y188L is selected by DOR and confers >50-fold reduced susceptibility to DOR (13,18,42).

Y188C is an uncommon nonpolymorphic mutation selected by NVP and EFV. It confers >30-fold reduced susceptibility to NVP (6,18,29,38). Its effect on EFV has varied between different studies. It does not appear to reduce susceptibility to ETR, RPV, or DOR.

Y188H is an uncommon nonpolymorphic mutation selected by NVP and EFV. It confers about 5-fold reduced susceptibility to these NRTIs (18,29,38). It does not confer reduced susceptibility to ETR, RPV, or DOR.

Y188F is a rare nonpolymorphic NNRTI-associated mutations that is usually present as part of an electrophoretic mixture. It likely represents a transitional mutation between Y and L or a partial revertant mutation.

G190A/S/E/Q
G190A is a nonpolymorphic mutation selected in persons receiving NVP and EFV (26,29,36). It confers >50-fold reduced NVP susceptibility and 5-10-fold reduced EFV susceptibility (68,74). It does not appear to be selected by RPV, ETR, or DOR or to reduce their antiviral activity (6,7,8,10,18).

G190S is a nonpolymorphic mutation selected in persons receiving NVP and EFV. It confers >50-fold decreased susceptibility to NVP and EFV (68,74). It does not appear to reduce susceptibility to ETR or RPV (6,18,69). It causes a variable reduction in DOR susceptibility ranging between 1.5 and 11-fold in clinical isolates (18). is the most commonNNRTI-resistance mutations in subtype A6 viruses from the countries of the former Soviet Union because the wildtype glycine in this lineage GGC requires just a single G-to-A mutation for G190S to develop (75).

G190E is an uncommon nonpolymorphic mutation associated with reduced replication capacity (74,76). It is reported in about 1% of persons receiving EFV (36,77) and in cell culture with ETR and RPV (1,2). It confers >100-fold reduced susceptibility to NVP and EFV and >10-fold reduced susceptibility to ETR and RPV (2,6,14,74). In the only isolate undergoing DOR susceptibility testing, it was associated with an 18-fold reduction in susceptibility (11,44). When present as a minority variant by NGS, it often occurs in the context of APOBEC-mediated G-to-A hypermutation and it is unlikely to be clinically relevant in this setting (78).

G190Q is an uncommon mutation occurring in about 0.5% of persons receiving EFV. It confers >100-fold reduced susceptibility to EFV and NVP. Its effects on ETR, RPV, and DOR susceptibility are not known. G190C/T/V are extremely rare nonpolymorphic NNRTI-selected mutations that confer high-level resistance to NVP and EFV (74). Their effects on ETR, RPV, and DOR susceptibility are not known.

H221Y
H221Y is a common nonpolymorphic accessory NNRTI-selected mutation that usually occurs in combination with Y181C (31,36,79). Alone it has minimal effects on NNRTI susceptibility (6). It is frequently selected in persons receiving RPV (3,4).

P225H
P225H is a nonpolymorphic accessory mutation selected in persons receiving EFV that usually occurs in combination with K103N (36,77,80). The combination of K103N and P225H confers >50-fold reduced susceptibility to NVP and EFV (6,68) and 5-10-fold reduced susceptibility to DOR (11,44).

F227L/C/I/V
F227L is a nonpolymorphic mutation selected in persons receiving NVP or EFV that usually occurs in combination with V106A (36,43). Alone F227L has little effect on NNRTI susceptibility. The combination of V106A and F227L confers high-level resistance to NVP and EFV (6,81) but not to ETR and RPV (8,24). The combination of V106A plus F227L has been selected in vitro by DOR (13,82) and is associated with >100-fold reduced DOR susceptibility (13,24).

F227C is a rare nonpolymorphic mutation that is usually selected in combination with other mutations often at position 106. It is selected in vitro by DOR (13,82) and in persons with HIVDR while receiving DOR (42,52). In combination with A98G, V106I, and V106M, it confers >100-fold reduced susceptibility to DOR (11,44). In combination with a variety of other mutations, it has generally been associated with high-level resistance to the remaining NNRTIs (2,6). F227I/V are extremely rare mutations that have been selected in vitro by DOR but have not been studied phenotypically (13).

M230L/I
M230L is an uncommon nonpolymorphic mutation selected in vitro by NVP, ETR, RPV, and DOR (1,2,14,82,83) and in vivo by NVP, EFV, and RPV (4,79). It is associated with reduced replication capacity (84). Alone, it causes variable reductions in NNRTI susceptibility – up to 5-fold for ETR and RPV and >10-fold for NVP, EFV, and DOR (2,14,24,84)

M230I is a rare NNRTI-resistance mutation selected in vitro by RPV (2,13). It appears to be associated with about 2-3-fold reduced susceptibility to ETR and RPV, 5-fold reduced susceptibility to EFV, and 10-fold reduced susceptibility to NVP (2,38,73). Its effect on DOR has not been studied. M230I is also commonly observed in sequences with APOBEC-mediated G-to-A hypermutation and in this setting it should not be considered clinically relevant (85).

Y232H
Y232H is a rare nonpolymorphic NNRTI-associated mutation that nearly always occurs in combination with other NNRTI-resistance mutations. It is selected in vitro by NVP and contributes to reduced NVP and EFV susceptibility (86,87).

L234I
L234I is an uncommon nonpolymorphic mutation selected in vitro by ETR (14) and DOR (13,82) and in persons receiving NVP and EFV. In combination with V106A, it is associated with >100-fold reduced susceptibility to DOR (13,24). Its effect on other NNRTIs is not known.

P236L
P236L is a rare nonpolymorphic mutation associated with high-level resistance to delavirdine (88,89). It confers about 4-fold reduced NVP susceptibility but does not appear to reduce susceptibility to other NNRTIs (2,18,38).

K238T/N
K238N/T are uncommon nonpolymorphic mutations selected by NVP and EFV usually in combination with K103N (36). Alone K103T reduces NVP susceptibility about 5 fold but does not appear to reduce EFV or ETR susceptibility (21,38). Alone, K238N appears to have minimal effect on NNRTI susceptibility.

Y318F
Y318F is a nonpolymorphic mutation selected in 2 of the 10 reported persons with VF and HIVDR while receiving a DOR-containing regimen (42) and in about 1% of persons receiving NVP and EFV (36). It has also been selected in vitro by ETR (14) and DOR (82). Alone, it was associated with a median 11-fold reduced DOR susceptibility in 3 clinical isolates (10) but otherwise it appears to have minimal if any effect on NVP, EFV, and ETR susceptibility (14,90).

N348I
N348I is a common nonpolymorphic accessory mutation in a region of RT that often doesn’t undergo sequencing. It appears to be selected by NVP, EFV, and NRTIs (91,92,93), Alone it reduces susceptibility to AZT, NVP, and EFV by about 2-fold (17,69,91,92). It has been proposed that N348I may reduce RNaseH cleavage thereby allowing more time for NNRTI dissociation and re-initiation of DNA synthesis (94,95).
"""

PI_STANFORD_TEXT = """L10F/I/V/R/Y
L10F is a common nonpolymorphic accessory mutation selected in persons receiving LPV, DRV, and ATV (11,12,13,14). It is associated with reduced in vitro susceptibility to LPV and DRV (15,16,17). In combination with other mutations, it may be associated with a reduced virological response to LPV and DOR (18,19,20).

L10I/V are polymorphic accessory PI-selected mutations that either reduce PI susceptibility or increase the replication of viruses containing PI-resistance mutations (21,22). L10R/Y are rare nonpolymorphic PI-selected mutations (23). Their effects on PI susceptibility have not been well studied.

V11I/L
V11I is a minimally polymorphic accessory mutation that is often selected in patients receiving DRV (24,25). In a weighted genotypic score, V11I was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26). V11L is a less common nonpolymorphic PI-selected mutation that may also be weakly associated with reduced DRV susceptibility (16,20).

K20R/I/M/T/V
K20T is a nonpolymorphic PI-selected accessory mutation (13,27,28,29). K20I is the consensus amino acid in subtype G and CRF02_AG viruses. However, in most other subtypes, it is a PI-selected mutation. K20M/V are uncommon relatively nonpolymorphic PI-selected mutations that have not been well studied (29). K20R is a highly polymorphic PI-selected accessory mutation that increases replication fitness in viruses with PI-resistance mutations (30).

L23I
L23I is an uncommon nonpolymorphic substrate cleft mutation selected primarily by NFV. It appears to only be associated with reduced NFV susceptibility (31).

L24I/F/M
L24I is a nonpolymorphic mutation selected primarily in persons receiving IDV and less often in persons receiving LPV (14,28,32). It is associated with reduced LPV and ATV susceptibility (15,16,33). L24F/M are rare PI-selected mutations.

D30N
D30N is a nonpolymorphic substrate-cleft mutation selected by NFV that reduces susceptibility to NFV but not to other PIs (15,16).

V32I
V32I is a nonpolymorphic substrate-cleft mutation selected in persons receiving LPV, ATV, and DRV (5,11,13,14,24,25,34,35). It reduces susceptibility to LPV, ATV, and DRV (15,16,26,33,36). In nearly two-thirds of sequences it occurs in combination with I47V/A (37,38). In combination with I47V it causes intermediate resistance to both LPV and DRV and lays the groundwork for higher levels of resistance with additional DRMs. In combination with I47A it causes high-level LPV resistance and intermediate DRV resistance. In a weighted genotypic score, V32I was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26).

L33F
L33F occurs in about 1.0% of subtype A, CRF01_AE, and CRF02_AG viruses from ART-naive persons but is otherwise nonpolymorphic. It is selected in persons receiving LPV, ATV, and DRV (5,11,13,14,24,25,34,35,39). When it occurs in combination with other PI-resistance mutations, it is associated with reduced susceptibility to LPV, ATV, and DRV (15,16,26). In a weighted genotypic score, L33F was one of 11 mutations associated with a reduced virological response to treatment with an DRV-containing regimen when combined with other DRV-resistance mutations (26).

L33I is a minimally polymorphic PI-selected mutation that does not appear to reduce PI susceptibility. L33V is a polymorphic mutation that is not selected by any of the PIs and does not reduce PI susceptibility.

M36I
M36I is the consensus amino acid in most non-B subtypes. It occurs in about 15% of PI-naive and 35% of PI-experienced individuals with subtype B viruses. In subtype B viruses, M36I increases the replication fitness of viruses with PI-resistance mutations (30,40).

K43T
K43T is a nonpolymorphic accessory mutation selected by multiple PIs including ATV and LPV (27,28,29). It was first recognized as a PI-resistance mutation because of its inclusion in a genotypic susceptibility score for TPV (41,42). It usually occurs with multiple other PI-resistance mutations making it difficult to discern its effect on PI susceptibility.

M46I/L/V
Among ART-naïve persons, M46I/L occur in 0.4% to 0.8% of CRF01_AE viruses and 0.1% to 0.3% of viruses belonging to other subtypes. M46I/L are selected by each of the PIs except SQV and DRV (13,14,28,35,39,43,44,45). M46I and M46L occur in about 20% and 10% of PI-treated patients, respectively. M46I usually occurs alone or in combination with V32I, I47V, L76V, I84V, and L90M whereas M46L usually occurs alone or in combination with I54V and V82A (38). M46I increases PR catalytic efficiency (46,47). M46I/L are associated with reduced susceptibility to ATV and LPV (15,16,33). M46V is an uncommon nonpolymorphic PI-selected mutation that has not been well studied.

I47V/A
I47V is a nonpolymorphic mutation selected by IDV, FPV, LPV and DRV (11,24,25,28,32,34,39,48,49). It is associated with reduced susceptibility to LPV and DRV but not to ATV (15,16,50). In a weighted genotypic score, I47V was one of 11 mutations associated with a reduced virological response to treatment with an DRV-containing regimen when combined with other DRV-resistance mutations (26).

I47A is a 2-base-pair nonpolymorphic mutation selected by LPV nearly always in combination with V32A (14,32,45,51,52). It confers high-level resistance to LPV and probably some degree of cross-resistance to DRV (15,16).

G48V/M/A/S/T/Q
G48V is a nonpolymorphic substrate-cleft mutation selected by SQV and less often by IDV and LPV (53,54). It confers intermediate resistance to ATV but has little if any effect on LPV susceptibility (15,16,55).

G48M is an uncommon 2-base-pair nonpolymorphic substrate-cleft mutation nearly always selected in viruses with multiple PI-resistance mutations. It has a resistance profile similar to G48V (15,16,50). G48A/S/T/Q/L are extremely rare nonpolymorphic PI-selected mutations nearly always selected in viruses with multiple PI-resistance mutations(29).

I50V/L
I50V is a nonpolymorphic substrate-cleft mutation selected by FPV, LPV and DRV (14,24,25,28,39). It reduces susceptibility to LPV and DRV (15,16,17,50). In a weighted genotypic score, I50V was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26).

I50L is a nonpolymorphic mutation selected by ATV (13,56,57,58,59). It confers high-level resistance to ATV and increases susceptibility to the remaining PIs (16,56,57,58,60).

F53L/Y
F53L is a nonpolymorphic accessory mutation selected primarily by SQV, IDV, ATV and LPV (13,28,39,44). In combination with other mutations, It is associated with reduced susceptibility to ATV and possibly LPV (15,16). F53Y is an uncommon nonpolymorphic accessory PI-selected mutation that has not been well studied (29).

I54V/A/S/T/L/M
I54V is a nonpolymorphic mutation selected primarily by IDV and LPV (28,32,34,39,43,54). It reduces susceptibility to each of the PIs except DRV (15,16,33).

I54M/L are nonpolymorphic mutations selected primarily by FPV, LPV and DRV (24,49,61). They are associated with reduced susceptibility to LPV, ATV, and DRV (15,16,17,50). In a weighted genotypic score, I54M/L were each one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26).

I54A/T/S are nonpolymorphic PI-selected mutations that occur exclusively in viruses with multiple PI-resistance mutations. They are associated with reduced susceptibility to each of the PIs except DRV (15,16).

Q58E
Q58E is a minimally polymorphic accessory mutation selected by each of the PIs except DRV. It was first recognized as a PI-resistance mutation because of its inclusion in a TPV genotypic susceptibility score (41,42). In combination with other PI-resistance mutations, it may contribute to low-level ATV resistance (15,16).

A71V/T/I/L
A71V/T are common polymorphic PI-selected accessory mutations that increase replication capacity and/or reduce PI susceptibility in viruses containing other PI-resistance mutations (30,62). A71I/L are nonpolymorphic accessory PI-selected mutations that occur in viruses containing multiple PI-resistance mutations (29).

G73S/T/C/A/D/V
G73S/T/C/A are common non-polymorphic accessory mutations selected primarily by SQV, ATV, IDV and NFV (28,39,43,53,56). They are associated with minimally reduced susceptibility to each of the PIs (15,16,17). In a weighted genotypic score, G73S was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26). However, it was replaced by T74P in an updated DRV genotypic score. G73V/D are rare nonpolymorphic PI-selected mutations (29).

T74P/S
T74P is a nonpolymorphic PI-selected accessory mutation that occurs primarily in viruses from persons who have received multiple PIs. It is associated with minimally reduced susceptibility to ATV and DRV (16). In a weighted genotypic score, it was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26). T74S is a PI-selected accessory mutation that is polymorphic in most non-B subtypes (63).

L76V
L76V is selected by IDV, LPV and DRV (14,25,35,64,65). It reduces susceptibility to LPV and DRV and increases susceptibility to ATV (15,16,17,66). In a weighted genotypic score, it was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26).

V82A/T/S/F/L/M/C
V82A is nonpolymorphic substrate-cleft mutation selected primarily by IDV and LPV (14,35,43,54,67). It is associated with reduced susceptibility to LPV and to a lesser extent ATV; it increases DRV susceptibility (15,16,33).

V82T/S are nonpolymorphic substrate-cleft mutations selected primarily by IDV and LPV (14,35,43,54,67). They are associated with reduced susceptibility to LPV and ATV, and increased susceptibility to DRV (15,16,33,50).

V82F is a nonpolymorphic substrate-cleft mutation selected primarily by IDV and LPV (33,39,43,45). It reduces LPV and DRV susceptibility (15,16,17).

V82L is a rare nonpolymorphic substrate-cleft mutation selected primarily by TPV (41). Its effect on other PIs is not well characterized.

V82M is an uncommon nonpolymorphic PI-selected substrate-cleft mutation. In most subtypes, it requires a 2-bp mutation and occurs in viruses from persons with multiple PI-resistance mutations. In subtype G, it usually requires a single bp mutation and has been reported it persons who have received IDV (68). It reduces susceptibility to IDV and possibly LPV (68).

V82C is an uncommon nonpolymorphic 2-base-pair PI-selected substrate-cleft mutation that occurs in viruses with multiple PI-resistance mutations. Its effect on PI susceptibility has not been well studied.

V82I is the consensus amino acid for position 82 in 90% of subtype G viruses and in 1% to 5% of the remaining subtypes. It does not appear to reduce PI susceptibility.

N83D/S
N83D is a nonpolymorphic PI-selected accessory mutation which possibly contributes to reduced ATV susceptibility in combination with other PI-resistance mutations (16). N83S is an extremely rare nonpolymorphic PI-selected mutation that has not been well studied.

I84V/A/C
I84V is a nonpolymorphic substrate-cleft mutation selected by each of the PIs (25,32,34,39,43,45,69,70,71). I84V reduces susceptibility to LPV, ATV, and DRV (15,16,17,33). In a weighted genotypic score, it was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26).

I84A/C are extremely rare non-polymorphic mutations selected primarily by SQV (72). I84A is associated with markedly reduced susceptibility to each of the PIs (15,72,73). I84C has a less-marked effect on PI susceptibility compared with I84A.

I85V
I85V is a nonpolymorphic PI-selected mutation (29). It has minimal, if any, effects on PI susceptibility (16).

N88D/S/T/G
N88S is a nonpolymorphic mutation selected by NFV, ATV, and occasionally IDV (6,13,59,74,75,76,77,78). It confers about 10-fold reduced ATV susceptibility (15,16,56,59).

N88D is a nonpolymorphic mutation selected by NFV, usually in combination with D30N (74,75,79). It is associated with cross-resistance to ATV (15,16).

N88G/T are extremely rare nonpolymorphic PI-selected mutations (29). Their effect on PI susceptibility appears to be minimal.

L89V/T
L89V is a nonpolymorphic accessory mutation weakly selected by each of the PIs. It appears to be minimally associated with reduced PI susceptibility (16). In a weighted genotypic score, it was one of 11 mutations associated with a reduced virological response to treatment with a DRV-containing regimen when combined with other DRV-resistance mutations (26).

L89T is a rare non-polymorphic PI-selected mutation selected primarily by ATV in the non-B subtypes for which the consensus mutation at position 89 is M (13) because in these subtypes, it is caused by a single C to T transition (ATG => ACG). Its effect on PI susceptibility has not been well studied.

L90M
L90M is selected primarily by SQV, NFV, IDV, and ATV (13,43,53,75) It contributes reduced susceptibility to ATV and to a lesser extent LPV (15,16)."""

rt_drm2polymorphic = {'M184V': NONPOLYMORPHIC, 'M184I': NONPOLYMORPHIC, 'K65R': NONPOLYMORPHIC,
                      'M41L': NONPOLYMORPHIC, 'D67N': NONPOLYMORPHIC, 'K70R': NONPOLYMORPHIC, 'L210W': NONPOLYMORPHIC,
                      'T215Y': NONPOLYMORPHIC, 'T215F': NONPOLYMORPHIC, 'K219Q': NONPOLYMORPHIC,
                      'K219E': NONPOLYMORPHIC,
                      'V118I': POLYMORPHIC, 'E40F': NONPOLYMORPHIC, 'E44A': NONPOLYMORPHIC, 'E44D': NONPOLYMORPHIC,
                      'K43Q': NONPOLYMORPHIC, 'K43N': NONPOLYMORPHIC, 'E203K': NONPOLYMORPHIC, 'H208Y': NONPOLYMORPHIC,
                      'D218E': NONPOLYMORPHIC, 'K223Q': NONPOLYMORPHIC, 'K223E': NONPOLYMORPHIC,
                      'L228H': NONPOLYMORPHIC, 'L228R': NONPOLYMORPHIC,
                      'A62V': NONPOLYMORPHIC, 'K70E': NONPOLYMORPHIC, 'K70G': NONPOLYMORPHIC, 'K70Q': NONPOLYMORPHIC,
                      'K70T': NONPOLYMORPHIC, 'K70N': NONPOLYMORPHIC,
                      'L74V': NONPOLYMORPHIC, 'L74I': NONPOLYMORPHIC, 'Y115F': NONPOLYMORPHIC, 'Q151M': NONPOLYMORPHIC,
                      'S68G': POLYMORPHIC,
                      'T69D': NONPOLYMORPHIC, 'T69N': NONPOLYMORPHIC, 'T69S': NONPOLYMORPHIC, 'T69G': NONPOLYMORPHIC,
                      'V75T': NONPOLYMORPHIC, 'V75M': NONPOLYMORPHIC, 'V75A': NONPOLYMORPHIC, 'V75S': NONPOLYMORPHIC,
                      'N348I': NONPOLYMORPHIC,
                      'T215S': NONPOLYMORPHIC, 'T215D': NONPOLYMORPHIC, 'K219R': NONPOLYMORPHIC,
                      'K219N': NONPOLYMORPHIC}

for line in NNRTI_STANFORD_TEXT.split('\n'):
    words = line.split(' ')
    if len(words) > 1:
        drm = words[0]
        drms = []
        if '/' in drm:
            start = drm.split('/')[0][:-1]
            for ending in drm[drm.find('/') - 1:]:
                if ending != '/':
                    drms.append(start + ending)
        else:
            drms = [drm]
        status = NONPOLYMORPHIC if NONPOLYMORPHIC in words else (
            POLYMORPHIC if POLYMORPHIC in words else (NONPOLYMORPHIC if 'selected' in line else 'unknown'))
        for drm in drms:
            rt_drm2polymorphic[drm] = status

drm2polymorphic = {}
for drm, pol in rt_drm2polymorphic.items():
    drm2polymorphic['RT_' + drm] = pol

for line in PI_STANFORD_TEXT.split('\n'):
    words = line.split(' ')
    if len(words) > 1:
        drm = words[0]
        drms = []
        if '/' in drm:
            start = drm.split('/')[0][:-1]
            for ending in drm[drm.find('/') - 1:]:
                if ending != '/':
                    drms.append(start + ending)
        else:
            drms = [drm]
        status = NONPOLYMORPHIC if NONPOLYMORPHIC in words else (
            POLYMORPHIC if POLYMORPHIC in words else (NONPOLYMORPHIC if 'selected' in line else 'unknown'))
        for drm in drms:
            drm2polymorphic['PR_' + drm] = status

drm2polymorphic['PR_M36I'] = POLYMORPHIC
drm2polymorphic['PR_V82I'] = POLYMORPHIC

drm2polymorphic['PR_M46I'] = NONPOLYMORPHIC
drm2polymorphic['PR_M46L'] = NONPOLYMORPHIC
drm2polymorphic['PR_M46V'] = NONPOLYMORPHIC


df = pd.DataFrame.from_dict(drm2polymorphic, orient='index')
df.columns = ['type']


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser(description='Annotate a tree with location data.')

    parser.add_argument('--tab', required=True, type=str,
                        help='Path to the output table.')
    params = parser.parse_args()

    df.to_csv(params.tab, sep='\t', index_label='DRM')
