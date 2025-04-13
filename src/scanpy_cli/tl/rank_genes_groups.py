import rich_click as click
import scanpy as sc
import sys
import pickle


@click.command(name="rank-genes-groups")
@click.option(
    "--groupby",
    type=str,
    required=True,
    help="The key of the observations grouping to consider.",
)
@click.option(
    "--groups",
    type=str,
    default="all",
    help="Subset of groups, e.g. ['g1', 'g2', 'g3'], to which comparison shall be restricted, or 'all' for all groups.",
)
@click.option(
    "--reference",
    type=str,
    default="rest",
    help="If 'rest', compare each group to the union of the rest of the group. If a group identifier, compare with respect to this group.",
)
@click.option(
    "--n-genes",
    type=int,
    help="The number of genes that appear in the returned tables (default: None).",
)
@click.option(
    "--use-raw",
    type=bool,
    default=None,
    help="Use raw attribute of adata if present (default: None).",
)
@click.option(
    "--layer",
    type=str,
    help="Key from adata.layers whose value will be used to perform tests on (default: None).",
)
@click.option(
    "--method",
    type=click.Choice(["logreg", "t-test", "wilcoxon", "t-test_overestim_var"]),
    help="The default method is 't-test'.",
)
@click.option(
    "--corr-method",
    type=click.Choice(["benjamini-hochberg", "bonferroni"]),
    default="benjamini-hochberg",
    help="Correction method for multiple testing (default: 'benjamini-hochberg').",
)
@click.option(
    "--mask-var",
    type=str,
    help="Select subset of genes to use in statistical tests (default: None).",
)
@click.option(
    "--rankby-abs",
    type=bool,
    default=False,
    help="Rank genes by the absolute value of the score, not by the score (default: False).",
)
@click.option(
    "--pts",
    type=bool,
    default=False,
    help="Compute the fraction of cells expressing the genes (default: False).",
)
@click.option(
    "--tie-correct",
    type=bool,
    default=False,
    help="Use tie correction for 'wilcoxon' scores (default: False).",
)
@click.option(
    "--input-file",
    "-i",
    required=True,
    help="Input h5ad file containing AnnData object.",
)
@click.option(
    "--output-file",
    "-o",
    required=True,
    help="Output h5ad file to save the processed AnnData object.",
)
@click.option(
    "--rank-genes-output",
    type=str,
    help="Optional path to save rank_genes_groups dictionary as a pickle file.",
)
def rank_genes_groups(
    groupby,
    groups,
    reference,
    n_genes,
    use_raw,
    layer,
    method,
    corr_method,
    mask_var,
    rankby_abs,
    pts,
    tie_correct,
    input_file,
    output_file,
    rank_genes_output,
):
    """Rank genes for characterizing groups.

    This function computes differential expression between groups using various statistical tests.
    It can be used to identify marker genes for different cell types or conditions.

    Results are stored in the AnnData object:
    - adata.uns['rank_genes_groups']: Dictionary containing the results
    - adata.uns['rank_genes_groups']['names']: Names of genes
    - adata.uns['rank_genes_groups']['scores']: Scores for each gene
    - adata.uns['rank_genes_groups']['pvals']: P-values for each gene
    - adata.uns['rank_genes_groups']['pvals_adj']: Adjusted p-values for each gene
    - adata.uns['rank_genes_groups']['pts']: Fraction of cells expressing the genes (if pts=True)
    - adata.uns['rank_genes_groups']['pts_rest']: Fraction of cells from the rest expressing the genes (if reference='rest')
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's rank_genes_groups function
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            groups=groups,
            reference=reference,
            n_genes=n_genes,
            use_raw=use_raw,
            layer=layer,
            method=method,
            corr_method=corr_method,
            mask_var=mask_var,
            rankby_abs=rankby_abs,
            pts=pts,
            tie_correct=tie_correct,
        )

        # Save the result
        adata.write(output_file)
        click.echo(f"Successfully ran rank_genes_groups and saved to {output_file}")

        # Save rank_genes_groups dictionary as pickle if specified
        if rank_genes_output:
            with open(rank_genes_output, "wb") as f:
                pickle.dump(adata.uns["rank_genes_groups"], f)
            click.echo(
                f"Successfully saved rank_genes_groups dictionary to {rank_genes_output}"
            )

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)
