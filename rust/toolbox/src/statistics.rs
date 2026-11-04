/// Shared Pearson calculation. Callers choose the measured slots and minimum
/// sample count; missing-slot policy is part of the scientific method.
pub(crate) fn pearson(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.len() < 2
        || left.len() != right.len()
        || left.iter().chain(right).any(|v| !v.is_finite())
    {
        return None;
    }
    let mean_left = left.iter().sum::<f64>() / left.len() as f64;
    let mean_right = right.iter().sum::<f64>() / right.len() as f64;
    let mut covariance = 0.0;
    let mut variance_left = 0.0;
    let mut variance_right = 0.0;
    for (left, right) in left.iter().zip(right) {
        let a = left - mean_left;
        let b = right - mean_right;
        covariance += a * b;
        variance_left += a * a;
        variance_right += b * b;
    }
    let correlation = covariance / (variance_left * variance_right).sqrt();
    (variance_left > 0.0 && variance_right > 0.0 && correlation.is_finite()).then_some(correlation)
}
