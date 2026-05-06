from __future__ import annotations


def clean_family_id(seq_id: str) -> str:
    """
    Normalize TE family identifiers.

    Examples
    --------
    Copia-2_Hm#LTR/Copia
        -> Copia-2_Hm

    CR1-1_Hmel_A#LINE/CR1
        -> CR1-1_Hmel_A
    """

    return seq_id.split("#", 1)[0]
