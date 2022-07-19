"""Analyze BioPAX models for enforcing reasonable prefix standards."""

from collections import Counter, defaultdict

import bioregistry
import click
import pybiopax
from pybiopax.biopax import BioPaxModel, Xref
from tabulate import tabulate


@click.command()
def main():
    """Run the example analysis of a BioPAX model."""
    model: BioPaxModel = pybiopax.model_from_reactome("418592")
    refs = Counter((ref.db, ref.id) for ref in model.get_objects_by_type(Xref))
    dd = defaultdict(dict)
    for (prefix, identifier), count in refs.items():
        dd[prefix][identifier] = count

    invalid = []
    unnorm = []
    for prefix in dd:
        norm_prefix = bioregistry.normalize_prefix(prefix)
        if norm_prefix is None:
            invalid.append(prefix)
        elif norm_prefix != prefix:
            unnorm.append((prefix, norm_prefix))
        else:
            pass

    print("Invalid\n")
    print(
        tabulate(
            Counter(invalid).most_common(),
            headers=["prefix", "norm_prefix", "count"],
            tablefmt="github",
        )
    )

    print("\nUnnormalized\n")
    print(
        tabulate(
            [(*k, v) for k, v in Counter(unnorm).most_common()],
            headers=["prefix", "count"],
            tablefmt="github",
        )
    )


if __name__ == "__main__":
    main()
