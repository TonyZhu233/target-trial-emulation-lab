"""Compatibility script for the original console-only demonstration."""

if __package__:
    from .core import analyze, build_cohort, make_synthetic_database
else:
    from core import analyze, build_cohort, make_synthetic_database


def main() -> None:
    db = make_synthetic_database()
    try:
        results = analyze(build_cohort(db))
    finally:
        db.close()
    for name, value in results.items():
        print(f"{name}: {value:.4f}")
    print("Full report: python -m src.run")


if __name__ == "__main__":
    main()
