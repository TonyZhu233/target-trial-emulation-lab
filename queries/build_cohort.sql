-- One row per eligible new user at time zero.
SELECT
    p.patient_id,
    p.age,
    p.comorbidity,
    CASE WHEN e.treatment = 'A' THEN 1 ELSE 0 END AS treated,
    CASE WHEN o.outcome_day BETWEEN 1 AND 180 THEN 1 ELSE 0 END AS outcome
FROM patients AS p
JOIN index_exposures AS e
    ON e.patient_id = p.patient_id
LEFT JOIN outcomes AS o
    ON o.patient_id = p.patient_id
WHERE p.age BETWEEN 50 AND 90
  AND e.prior_enrollment_days >= 365
  AND e.prior_class_use = 0;

