# Arkworks 0.5.0 兼容性修复

本文档总结了为使 poly-commit 库与 arkworks 0.5.0 兼容而进行的修复。

## 主要问题

1. **Sync trait 缺失**: arkworks 0.5.0 要求某些参数类型实现 Sync trait
2. **Polynomial trait 缺失**: MultilinearExtension 类型需要同时实现 Polynomial trait

## 修复的文件

### 1. linear_codes/brakedown.rs
- 在所有 `PCUniversalParams`, `PCCommitterKey`, `PCVerifierKey` 实现中添加了 Sync 约束：
  - `<<C as Config>::LeafHash as CRHScheme>::Parameters: Sync`
  - `<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Sync`
  - `<H as CRHScheme>::Parameters: Sync`

### 2. linear_codes/ligero.rs
- 在所有 `PCUniversalParams`, `PCCommitterKey`, `PCVerifierKey` 实现中添加了相同的 Sync 约束

### 3. linear_codes/multilinear_brakedown/mod.rs
- 在 `LinearEncode` 实现中添加了：
  - `P: MultilinearExtension<F> + ark_poly::Polynomial<F>`
  - Sync 约束
- 更新了类型引用：`<P as ark_poly::Polynomial<F>>::Point`

### 4. linear_codes/multilinear_ligero/mod.rs
- 在 `LinearEncode` 实现中添加了相同的约束和类型引用更新

### 5. linear_codes/univariate_ligero/mod.rs
- 在 `LinearEncode` 实现中添加了 Sync 约束

### 6. hyrax/mod.rs
- 添加了 `Polynomial` trait 导入：`crate::Polynomial`
- 在 `PolynomialCommitment` 实现中添加了：
  - `P: MultilinearExtension<G::ScalarField> + ark_poly::Polynomial<G::ScalarField>`

## 修复的核心原理

1. **Sync 约束**: arkworks 0.5.0 中的序列化要求参数类型是线程安全的，因此必须实现 Sync trait。

2. **Polynomial trait**: 在 arkworks 0.5.0 中，`MultilinearExtension` 不再自动提供 `Polynomial` trait 的功能，需要显式地约束类型同时实现两个 trait。

3. **完整路径**: 使用 `ark_poly::Polynomial<F>` 而不是简单的 `Polynomial<F>` 来避免命名冲突。

## 验证

修复完成后，可以通过以下命令验证：

```bash
cargo check
cargo test
cargo bench
```

这些修复确保了 poly-commit 库与 arkworks 0.5.0 的完全兼容性。
