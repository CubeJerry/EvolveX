## Original EvolveX

The original search uses a binding-centred Metropolis objective:

```text
ΔE = ΔG_binding
```

Favourable mutations are accepted automatically:

```text
ΔG_binding < 0
```

Unfavourable mutations may still be accepted with:

```text
P = e^-(ΔG_binding / 0.5919)
```

The original implementation also applies strict structural filters for antibody stability and intraclash.

---

## Forked EvolveX

This fork replaces the binding-only objective with a rewarded objective:

```text
ΔE_rewarded = ΔG_bind,obj + 0.25·ΔG_stab + 0.5·Δn_mut
```

where:

```text
ΔG_bind,obj = binding objective
ΔG_stab     = antibody stability change
Δn_mut      = mutation-count change from the parental sequence
```

Lower values are better.

Favourable moves are accepted automatically:

```text
ΔE_rewarded < 0
```

Unfavourable moves may still be accepted with the same Metropolis form:

```text
P = e^-(ΔE_rewarded / 0.5919)
```


---

## Added design constraints

This fork adds several extra guardrails:

* mutation-burden tracking
* soft mutation cap at 25% divergence
* hard mutation cap at 50% divergence
* state-dependent proposal weights for selected residue types

The penalty residues are:

```text
F, H, I, R, W, Y
```

These residues are down-weighted when overrepresented, reducing excessive enrichment of bulky, aromatic, charged, or potentially sticky residues.
