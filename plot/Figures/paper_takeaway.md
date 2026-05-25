# Paper Takeaway — Bulk Omics Experiment (Fig.4)

1. **guided-PLSはbulk RNA-seqの種間diagonal integrationでも有効**。共通遺伝子空間を持たないhost/parasiteデータに対して、metadata-guided alignmentにより解釈可能な潜在因子を抽出できる。

2. **因子特異性(factor specificity)でguided-PLSが最高性能**。各次元が特定の実験条件に1対1対応する度合いが0.77–0.85と、他手法(0.50–0.81)を上回る。特にRBH-PLSは全次元が同一条件に固着(n_distinct=1)。

3. **種間concordanceがguided-PLSのみ100%**。AtとPjで同じ次元が同じ実験条件に対応する割合が4/4(100%)。guided-PCAは同じラベルを使うにもかかわらず14–55%に留まり、種間統合(PLS構造)とラベルガイドの両方が必要であることを示す。

4. **空間オミクス・single-cellに続くrobustness実証**。bulk RNA-seqという異なるモダリティ・スケールでもguided-PLSの優位性が再現され、手法のgeneralityを裏付ける。
