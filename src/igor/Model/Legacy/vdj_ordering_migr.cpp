#include <igor/Model/Legacy/vdj_ordering_migr.h>

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace igor::model::legacy {
namespace migration {

void require_vdj_or_vj_ordering(const igor::core::SeqTypeRegistry &registry, const std::string &reader)
{
  static const std::vector<std::string> vdj = {"V_gene_seq", "VD_ins_seq", "D_gene_seq", "DJ_ins_seq", "J_gene_seq"};
  static const std::vector<std::string> vj = {"V_gene_seq", "VJ_ins_seq", "J_gene_seq"};
  const std::vector<std::string> &ordering = registry.get_ordered_types();
  if (ordering == vdj or ordering == vj) {
    return;
  }

  //Name where the ordering leaves the shape it follows longest.
  const auto common_prefix = [&ordering](const std::vector<std::string> &shape) {
    return static_cast<std::size_t>(
        std::mismatch(ordering.begin(), ordering.end(), shape.begin(), shape.end()).first - ordering.begin());
  };
  const std::vector<std::string> &shape = common_prefix(vdj) >= common_prefix(vj) ? vdj : vj;
  const std::size_t at = common_prefix(shape);

  std::string listed;
  for (const std::string &name : ordering) {
    listed += (listed.empty() ? "" : ", ") + name;
  }
  throw std::invalid_argument(reader + " reads a V-D-J or a V-J model only, by the legacy segment names; this model's "
                              + "ordering is " + (listed.empty() ? std::string("empty") : listed) + ", "
                              + (at < ordering.size() ? "where it cannot place " + ordering[at]
                                                      : "which lacks " + shape[at]));
}

} // namespace migration
} // namespace igor::model::legacy
